/**************************************************************************
 * Copyright (c) 2019- Guillermo A. Perez
 * Copyright (c) 2020- Tom van Dijk
 * 
 * This file is a modified version of HOA2PG of HOATOOLS.
 * 
 * HOATOOLS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * HOATOOLS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with HOATOOLS. If not, see <http://www.gnu.org/licenses/>.
 * 
 * Guillermo A. Perez
 * University of Antwerp
 * guillermoalberto.perez@uantwerpen.be
 *
 * Tom van Dijk
 * University of Twente
 * t.vandijk@utwente.nl
 *************************************************************************/

#include <cassert> // for assert
#include <cstring> // for strcmp
#include <iostream>
#include <sys/time.h> // for gettimeofday

#include <game.hpp>
#include <oink.hpp>
#include <solvers.hpp>
#include <tools/cxxopts.hpp>
#include <sylvan.h>
#include <map>
#include <set>

extern "C" {
#include "simplehoa.h"
}

using namespace sylvan;

static double
wctime()
{
    struct timeval time;
    gettimeofday(&time, NULL);
    return time.tv_sec + 1E-6 * time.tv_usec;
}

/**
 * evalLabel converts a label on a transition to a MTBDD encoding the label
 */
#define evalLabel(a,b,c) CALL(evalLabel, a, b, c)
TASK_3(MTBDD, evalLabel, BTree*, label, AliasList*, aliases, uint32_t*, variables)
{
    MTBDD left = mtbdd_false, right = mtbdd_false, result = mtbdd_false;
    mtbdd_refs_pushptr(&left);
    mtbdd_refs_pushptr(&right);
    switch (label->type) {
        case NT_BOOL:
            result = label->id ? mtbdd_true : mtbdd_false;
            break;
        case NT_AND:
            left = evalLabel(label->left, aliases, variables);
            right = evalLabel(label->right, aliases, variables);
            result = sylvan_and(left, right);
            break;
        case NT_OR:
            left = evalLabel(label->left, aliases, variables);
            right = evalLabel(label->right, aliases, variables);
            result = sylvan_or(left, right);
            break;
        case NT_NOT:
            left = evalLabel(label->left, aliases, variables);
            result = sylvan_not(left);
            break;
        case NT_AP:
            result = mtbdd_ithvar(variables[label->id]);
            break;
        case NT_ALIAS:
            for (AliasList* a = aliases; a != NULL; a = a->next) {
                if (strcmp(a->alias, label->alias) == 0) {
                    return evalLabel(a->labelExpr, aliases, variables);
                }
            }
            break;
        default:
            assert(false);  // all cases should be covered above
    }
    mtbdd_refs_popptr(2);
    return result;
}




/**
 * This helper function ensures that the priority p is adjusted to
 * ensure we have a "max, even" parity game, as this is what Oink expects.
 */
static inline int
adjustPriority(int p, bool maxPriority, bool controllerIsOdd, int noPriorities)
{
    // To deal with max vs min, we subtract from noPriorities if
    // originally it was min (for this we need it to be even!)
    if (!maxPriority) {
        int evenMax = 2*((noPriorities+1)/2);
        p = evenMax - p;  // flip from min to max
    }
    // We reserve priority 0, so do not use it.
    p += 2;
    // If "controllerParity" is 1 (odd), automatically adjust the priority
    if (controllerIsOdd) p--;
    return p;
}




/**
 * Collect all intermediary MTBDD roots that become the vertices
 * where player 0 chooses the response in controllable APs
 */
void
collect_inter(MTBDD trans, int uap_count, std::set<MTBDD> &res)
{
    if (mtbdd_isleaf(trans)) {
        res.insert(trans);
    } else {
        // node
        uint32_t var = mtbdd_getvar(trans);
        if (var < (unsigned)uap_count) {
            // uncontrollable ap
            collect_inter(mtbdd_gethigh(trans), uap_count, res);
            collect_inter(mtbdd_getlow(trans), uap_count, res);
        } else {
            // controllable ap
            res.insert(trans);
        }
    }
}


/**
 * Given some intermediary MTBDD root, collect all the MTBDD
 * leafs, representing the target vertices of the full transition
 * and the transition priority encoded in a single 64-bit value
 */
void
collect_targets(MTBDD trans, std::set<uint64_t> &res)
{
    if (mtbdd_isleaf(trans)) {
        res.insert(mtbdd_getint64(trans));
    } else {
        collect_targets(mtbdd_gethigh(trans), res);
        collect_targets(mtbdd_getlow(trans), res);
    }
}



/**
 * Construct and solve the game explicitly
 */
pg::Game*
constructGame(HoaData *data, bool isMaxParity, bool controllerIsOdd)
{
    // Set which APs are controllable in the bitset controllable
    pg::bitset controllable(data->noAPs);
    for (IntList* c = data->cntAPs; c != NULL; c = c->next) controllable[c->i] = 1;

    // Count the number of controllable/uncontrollable APs
    const int cap_count = controllable.count();
    const int uap_count = data->noAPs - cap_count;

    // Initialize the BDD variable indices
    // Variable order uncontrollable < controllable
    uint32_t variables[data->noAPs];
    uint32_t uidx = 0, oidx = data->noAPs - controllable.count();
    for (int i=0; i<data->noAPs; i++) {
        if (controllable[i]) variables[i] = oidx++;
        else variables[i] = uidx++;
    }

    // Now initialize a new parity game
    pg::Game *game = new pg::Game(data->noStates * 10); // start with 10 the number of states
    // notice that the number of vertices automatically grows when needed anyway

    int nextIndex = data->noStates; // index of the next vertex to make

    std::vector<int> succ_state;  // for current state, the successors
    std::vector<int> succ_inter;  // for current intermediate state, the successors

    MTBDD trans_bdd = mtbdd_false;
    mtbdd_refs_pushptr(&trans_bdd);

    std::set<MTBDD> inter_bdds;
    std::set<uint64_t> targets;

    // these variables are used deeper in the for loops, however we can
    // push them to mtbdd refs here and avoid unnecessary pushing and popping
    MTBDD lblbdd = mtbdd_false;
    MTBDD leaf = mtbdd_false;
    mtbdd_refs_pushptr(&lblbdd);
    mtbdd_refs_pushptr(&leaf);

    LACE_ME;

    // Loop over every state
    for (StateList* state = data->states; state != NULL; state = state->next) {
        trans_bdd = mtbdd_false;

        for (TransList* trans = state->transitions; trans != NULL; trans = trans->next) {
            // there should be a single successor per transition
            assert(trans->successors != NULL && trans->successors->next == NULL);
            // there should be a label at state or transition level
            BTree* label;
            if (state->label != NULL) label = state->label;
            else label = trans->label;
            assert(label != NULL);
            // there should be a priority at state or transition level
            int priority = 0;
            if (state->accSig == NULL) {
                // there should be exactly one acceptance set!
                IntList* acc = trans->accSig;
                assert(acc != NULL && acc->next == NULL);
                // adjust priority
                priority = adjustPriority(acc->i, isMaxParity, controllerIsOdd, data->noAccSets);
            }            
            // tricky tricky
            lblbdd = evalLabel(label, data->aliases, variables);
            leaf = mtbdd_int64(((uint64_t)priority << 32) | (uint64_t)(trans->successors->i));
            trans_bdd = mtbdd_ite(lblbdd, leaf, trans_bdd);
            lblbdd = leaf = mtbdd_false;
        }

        // At this point, we have the transitions from the state all in a neat
        // single BDD. Time to generate the split game fragment from the current
        // state.

        collect_inter(trans_bdd, uap_count, inter_bdds);
        for (MTBDD inter_bdd : inter_bdds) {
            collect_targets(inter_bdd, targets);
            for (uint64_t lval : targets) {
                int priority = (int)(lval >> 32);
                int target = (int)(lval & 0xffffffff);

                if (priority != 0) {
                    int vfin = nextIndex++;
                    game->init_vertex(vfin, priority, 0);
                    game->e_start(vfin);
                    game->e_add(vfin, target);
                    game->e_finish();
                    succ_inter.push_back(vfin);
                } else {
                    succ_inter.push_back(target);
                }
            }

            int vinter = nextIndex++;
            game->init_vertex(vinter, 0, 0);
            game->e_start(vinter);
            for (int to : succ_inter) game->e_add(vinter, to);
            game->e_finish();
            succ_state.push_back(vinter);

            succ_inter.clear();
            targets.clear();
        }

        // there should be a priority at state or transition level
        int priority;
        if (state->accSig != NULL) priority = adjustPriority(state->accSig->i, isMaxParity, controllerIsOdd, data->noAccSets);
        else priority = 0;

        game->init_vertex(state->id, priority, 1, state->name ? state->name : "");
        game->e_start(state->id);
        for (int to : succ_state) game->e_add(state->id, to);
        game->e_finish();

        succ_state.clear();
        inter_bdds.clear();
    }

    // tell Oink we're done adding stuff, resize game to final size
    game->v_resize(nextIndex);

    mtbdd_refs_popptr(3);

    return game;
}


cxxopts::ParseResult
handleOptions(int &argc, char**& argv)
{
    try {
        cxxopts::Options opts(argv[0], "HOA synthesis using Sylvan and Oink");
        opts.custom_help("[OPTIONS...] [FILE]");
        opts.add_options()
            ("help", "Print help")
            ("symbolic", "Generate a symbolic parity game")
            ("print-game", "Just print the parity game")
            ("v,verbose", "Be verbose")
            ;
        opts.add_options("Explicit solvers")
            ("solvers", "List available solvers")
            ;

        // Add solvers
        pg::Solvers solvers;
        for (unsigned id=0; id<solvers.count(); id++) {
            opts.add_options("Explicit solvers")(solvers.label(id), solvers.desc(id));
        }

        // Parse command line
        auto options = opts.parse(argc, argv);

        if (options.count("help")) {
            std::cout << opts.help() << std::endl;
            exit(0);
        }

        if (options.count("solvers")) {
            solvers.list(std::cout);
            exit(0);
        }

        return options;
    } catch (const cxxopts::OptionException& e) {
        std::cout << "error parsing options: " << e.what() << std::endl;
        exit(0);
    }
}


/**
 * The main function
 */
int
main(int argc, char* argv[])
{
    auto options = handleOptions(argc, argv);

    bool verbose = options["verbose"].count() > 0;

    // First initialize the HOA data structure
    HoaData* data = (HoaData*)malloc(sizeof(HoaData));
    defaultsHoa(data);

    if (argc == 1) {
        int ret = parseHoa(stdin, data);
        if (ret != 0) return ret;
    } else {
        FILE* f = fopen(argv[1], "r");
        if (f == NULL) {
            std::cout << "file not found: " << argv[1] << std::endl;
            return 0;
        }
        int ret = parseHoa(f, data);
        fclose(f);
        if (ret != 0) return ret;
    }

    if (verbose) std::cerr << "finished reading file." << std::endl;

    // First check if the automaton is a parity automaton
    bool isMaxParity = true;
    short controllerParity = 0;
    int ret = isParityGFG(data, &isMaxParity, &controllerParity);
    if (ret != 0) return ret;
    bool controllerIsOdd = controllerParity != 0;

    // Initialize Lace
    lace_init(1, 0); // initialize Lace, but sequentially
    lace_startup(0, 0, 0); // no thread spawning

    // And initialize Sylvan
    sylvan_set_limits(128LL << 20, 1, 16); // should be enough (128 megabytes)
    sylvan_init_package();
    sylvan_init_mtbdd();

    bool explicit_solver = true;
    bool write_pg = options["print-game"].count() > 0;

    if (explicit_solver) {
        // Remember the start vertex
        int vstart = data->start->i;

        // Construct the game
        pg::Game *game = constructGame(data, isMaxParity, controllerIsOdd);
        game->set_label(vstart, "initial");

        // free HOA allocated data structure
        deleteHoa(data);
        if (verbose) std::cerr << "finished constructing game." << std::endl;

        // We don't need Sylvan anymore at this point
        sylvan_quit();

        if (write_pg) {
            // in case we want to write the file to PGsolver file format...
            game->write_pgsolver(std::cout);
            // std::cerr << "initial vertex: " << vstart << std::endl;
            exit(0);
        }

        // we sort now, so we can track the initial state
        int *mapping = new int[game->vertexcount()];
        game->sort(mapping);
        for (int i=0; i<game->vertexcount(); i++) {
            if (mapping[i] == vstart) {
                vstart = i;
                break;
            }
        }
        delete[] mapping;

        // OK, fire up the engine
       
        std::stringstream log;

        std::string solver = "fpi";
        pg::Solvers solvers;
        for (unsigned id=0; id<solvers.count(); id++) {
            if (options.count(solvers.label(id))) solver = solvers.label(id);
        }
 
        pg::Oink engine(*game, verbose ? std::cerr : log);
        engine.setTrace(0);
        engine.setRenumber();
        engine.setSolver(solver);
        engine.setWorkers(-1);

        // and run the solver
        double begin = wctime();
        engine.run();
        double end = wctime();

        // report how long it all took
        if (verbose) std::cerr << "total solving time: " << std::fixed << (end-begin) << " sec." << std::endl;

        // finally, check if the initial vertex is won by controller or environment
        if (game->winner[vstart] == 0) {
            std::cout << "REALIZABLE";
            exit(10);
        } else {
            std::cout << "UNREALIZABLE";
            exit(20);
        }
    }
}

