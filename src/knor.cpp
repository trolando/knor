/**************************************************************************
 * Copyright Tom van Dijk
 *************************************************************************/

#include <cassert> // for assert
#include <cstring> // for strcmp
#include <iostream>
#include <map>
#include <set>
#include <sys/time.h> // for gettimeofday
#include <deque>

#include <oink/game.hpp>
#include <oink/oink.hpp>
#include <oink/solvers.hpp>
#include <cxxopts.hpp>
#include <knor.hpp>
#include <symgame.hpp>
#include <bisim.hpp>
#include <aigmaker.hpp>

extern "C" {
    #include <sylvan.h>
    #include "simplehoa.h"
    #include "aiger.h"
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
 * Convert a transition label (Btree) to a BDD encoding the label
 * a label is essentially a boolean combination of atomic propositions and aliases
 */
TASK_IMPL_3(MTBDD, evalLabel, BTree*, label, HoaData*, data, uint32_t*, variables)
{
    MTBDD left = mtbdd_false;
    MTBDD right = mtbdd_false;
    MTBDD result = mtbdd_false;
    switch (label->type) {
        case NT_BOOL:
            return label->id ? mtbdd_true : mtbdd_false;
        case NT_AND:
            mtbdd_refs_pushptr(&left);
            mtbdd_refs_pushptr(&right);
            left = CALL(evalLabel, label->left, data, variables);
            right = CALL(evalLabel, label->right, data, variables);
            result = sylvan_and(left, right);
            mtbdd_refs_popptr(2);
            return result;
        case NT_OR:
            mtbdd_refs_pushptr(&left);
            mtbdd_refs_pushptr(&right);
            left = CALL(evalLabel, label->left, data, variables);
            right = CALL(evalLabel, label->right, data, variables);
            result = sylvan_or(left, right);
            mtbdd_refs_popptr(2);
            return result;
        case NT_NOT:
            mtbdd_refs_pushptr(&left);
            left = CALL(evalLabel, label->left, data, variables);
            result = sylvan_not(left);
            mtbdd_refs_popptr(1);
            return result;
        case NT_AP:
            return mtbdd_ithvar(variables[label->id]);
        case NT_ALIAS:
            // apply the alias
            for (int i=0; i<data->noAliases; i++) {
                Alias *a = data->aliases+i;
                if (strcmp(a->alias, label->alias) == 0) {
                    return CALL(evalLabel, a->labelExpr, data, variables);
                }
            }
            return mtbdd_invalid;
        default:
            assert(false);  // all cases should be covered above
            return mtbdd_invalid;
    }
}


/**
 * Given a label and a valuation of some of the atomic propositions,
 * we determine whether the label is true (1), false (-1), or its
 * value is unknown (0). The valuation is expected as an unsigned
 * integer whose i-th bit is 1 iff the i-th AP in apIds is set to 1
 * Returns -2 if there is a problem.
 */
static int
evalLabelNaive(BTree* label, Alias* aliases, int numAliases, int numAPs, int* apIds, uint64_t value) {
    assert(label != NULL);
    int left;
    int right;
    uint64_t mask;
    switch (label->type) {
        case NT_BOOL:
            return label->id ? 1 : -1;  // 0 becomes -1 like this
        case NT_AND:
            left = evalLabelNaive(label->left, aliases, numAliases, numAPs, apIds, value);
            right = evalLabelNaive(label->right, aliases, numAliases, numAPs, apIds, value);
            if (left == -1 || right == -1) return -1;
            if (left == 0 || right == 0) return 0;
            // otherwise
            return 1;
        case NT_OR:
            left = evalLabelNaive(label->left, aliases, numAliases, numAPs, apIds, value);
            right = evalLabelNaive(label->right, aliases, numAliases, numAPs, apIds, value);
            if (left == 1 || right == 1) return 1;
            if (left == 0 || right == 0) return 0;
            // otherwise
            return -1;
        case NT_NOT:
            return -1 * evalLabelNaive(label->left, aliases, numAliases, numAPs, apIds, value);
        case NT_AP:
            mask = 1;
            for (int i = 0; i < numAPs; i++) {
                if (label->id == apIds[i]) {
                    return ((mask & value) == mask) ? 1 : -1;
                }
                mask = mask << 1;
            }
            return 0;
        case NT_ALIAS:
            for (int i=0; i<numAliases; i++) {
                Alias* a = aliases+i;
                if (strcmp(a->alias, label->alias) == 0) {
                    return evalLabelNaive(a->labelExpr, aliases, numAliases, numAPs, apIds, value);
                }
            }
            break;
        default:
            assert(false);  // all cases should be covered above
    }
    return -2;
}



/**
 * This helper function ensures that the priority p is adjusted to
 * ensure we have a "max, even" parity game, as this is what Oink expects.
 * @param controllerIsOdd if the controller is the odd player
 * @param noPriorities how many priorities are in the game
 * @param maxPriority if the game is a max game
 */
int
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
 * where player 0 chooses the response in controllable APs.
 *
 * The input <trans> is the [partial] transition, essentially we just skip
 * the first <uap_count> variables, assuming that those are the uncontrollable APs...
 *
 * This is also why u_ap before c_ap variables in the BDD!
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
 * Encode priority i.e. all states via priority <priority>
 */
MTBDD
encode_prio(int priority, int priobits)
{
    MTBDD cube = mtbdd_true;
    for (int i=0; i<priobits; i++) {
        if (priority & 1) cube = mtbdd_makenode(priobits-i-1, mtbdd_false, cube);
        else cube = mtbdd_makenode(priobits-i-1, cube, mtbdd_false);
        priority >>= 1;
    }
    return cube;
}


/**
 * Encode a state as a BDD, using statebits 0..<statebits>, offsetted by <offset>+<priobits>
 * High-significant bits come before low-significant bits in the BDD
 */
MTBDD
encode_state(uint32_t state, const int statebits, const int s_first_var)
{
    // create a cube
    MTBDD cube = mtbdd_true;
    for (int i=0; i<statebits; i++) {
        const int bit = s_first_var+statebits-i-1;
        if (state & 1) cube = mtbdd_makenode(bit, mtbdd_false, cube);
        else cube = mtbdd_makenode(bit, cube, mtbdd_false);
        state >>= 1;
    }
    return cube;
}


/**
* Encode a priostate as a BDD, with priobits before statebits
* High-significant bits come before low-significant bits in the BDD
*/
MTBDD
encode_priostate(uint32_t state, uint32_t priority, const int statebits, const int priobits, const int s_first_var, const int p_first_var)
{
    // create a cube
    MTBDD cube = mtbdd_true;
    for (int i=0; i<statebits; i++) {
        if (state & 1) cube = mtbdd_makenode(s_first_var+statebits-i-1, mtbdd_false, cube);
        else cube = mtbdd_makenode(s_first_var+statebits-i-1, cube, mtbdd_false);
        state >>= 1;
    }
    for (int i=0; i<priobits; i++) {
        if (priority & 1) cube = mtbdd_makenode(p_first_var+priobits-i-1, mtbdd_false, cube);
        else cube = mtbdd_makenode(p_first_var+priobits-i-1, cube, mtbdd_false);
        priority >>= 1;
    }
    return cube;
}


/**
 * Given some intermediary MTBDD root, collect all the MTBDD leafs.
 * It decodes the leaf <priority, state> and re-encodes the leaf as a BDD such that
 * the first <priobits> BDD variables encode the priority and
 * the next <statebits> BDD variables encode the target state
 */
MTBDD
collect_targets(MTBDD trans, std::set<uint64_t> &res, const int statebits, const int priobits)
{
    if (mtbdd_isleaf(trans)) {
        uint64_t leaf = mtbdd_getint64(trans);
        res.insert(leaf);

        uint32_t priority = (uint32_t)(leaf>>32);
        uint32_t state = (uint32_t)(leaf & 0xffffffff);

        return encode_priostate(state, priority, statebits, priobits, priobits, 0);
    } else {
        MTBDD left = mtbdd_false;
        MTBDD right = mtbdd_false;
        mtbdd_refs_pushptr(&left);
        mtbdd_refs_pushptr(&right);

        left = collect_targets(mtbdd_getlow(trans), res, statebits, priobits);
        right = collect_targets(mtbdd_gethigh(trans), res, statebits, priobits);
        MTBDD result = sylvan_or(left, right);

        mtbdd_refs_popptr(2);
        return result;
    }
}


/**
 * Construct and solve the game explicitly
 */
TASK_3(pg::Game*, constructGameNaive, HoaData*, data, bool, isMaxParity, bool, controllerIsOdd)
{
    // Set which APs are controllable in the bitset controllable
    pg::bitset controllable(data->noAPs);
    for (int i=0; i<data->noCntAPs; i++) {
        controllable[data->cntAPs[i]] = 1;
    }

    // Count the number of controllable/uncontrollable APs
    const int cap_count = controllable.count();
    const int uap_count = data->noAPs - cap_count;
    const uint64_t numValuations = (1ULL << uap_count);

    int ucntAPs[uap_count];
    int uidx = 0;
    for (int i = 0; i < data->noAPs; i++) {
        if (!controllable[i]) ucntAPs[uidx++] = i;
    }

    // Now initialize a new parity game
    pg::Game *game = new pg::Game(data->noStates * 10); // start with 10 the number of states
    // notice that the number of vertices automatically grows when needed anyway

    int nextIndex = data->noStates; // index of the next vertex to make

    std::vector<int> succ_state;  // for current state, the successors
    std::vector<int> succ_inter;  // for current intermediate state, the successors

    // Loop over every state
    for (int i=0; i<data->noStates; i++) {
        auto state = data->states+i;
        for (uint64_t value = 0; value < numValuations; value++) {
            // for every valuation to the uncontrollable APs, we make an intermediate vertex
            for (int j=0; j<state->noTrans; j++) {
                auto trans = state->transitions+j;
                // there should be a single successor per transition
                assert(trans->noSucc == 1);
                // there should be a label at state or transition level
                BTree* label;
                if (state->label != NULL) label = state->label;
                else label = trans->label;
                assert(label != NULL);
                // we add a vertex + edges if the transition is compatible with the
                // valuation we are currently considering
                int evald = evalLabelNaive(label, data->aliases, data->noAliases, uap_count, ucntAPs, value);
                if (evald == -1) continue; // not compatible
                // there should be a priority at state or transition level
                if (state->accSig == NULL) {
                    // there should be exactly one acceptance set!
                    assert(trans->noAccSig == 1);
                    // adjust priority
                    int priority = adjustPriority(trans->accSig[0], isMaxParity, controllerIsOdd, data->noAccSets);

                    int vfin = nextIndex++;
                    game->init_vertex(vfin, priority, 0);
                    game->e_start(vfin);
                    game->e_add(vfin, trans->successors[0]);
                    game->e_finish();
                    succ_inter.push_back(vfin);
                } else {
                    succ_inter.push_back(trans->successors[0]);
                }
            }

            int vinter = nextIndex++;
            succ_state.push_back(vinter);
            game->init_vertex(vinter, 0, 0);
            game->e_start(vinter);
            for (int to : succ_inter) game->e_add(vinter, to);
            game->e_finish();
            succ_inter.clear();
        }

        int priority;
        if (state->accSig != NULL) {
            priority = adjustPriority(state->accSig[0], isMaxParity, controllerIsOdd, data->noAccSets);
        } else {
            priority = 0;
        }

        game->init_vertex(state->id, priority, 1, state->name ? state->name : std::to_string(state->id));
        game->e_start(state->id);
        for (int to : succ_state) game->e_add(state->id, to);
        game->e_finish();
        succ_state.clear();
    }

    // tell Oink we're done adding stuff, resize game to final size
    game->v_resize(nextIndex);

    return game;
}


/**
 * Construct and solve the game explicitly
 */
TASK_3(pg::Game*, constructGame, HoaData *, data, bool, isMaxParity, bool, controllerIsOdd)
{
    // Set which APs are controllable in the bitset controllable
    pg::bitset controllable(data->noAPs);
    for (int i=0; i<data->noCntAPs; i++) {
        controllable[data->cntAPs[i]] = 1;
    }

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
    // The number of vertices automatically grows when needed, but 10xstates is a good start

    int nextIndex = data->noStates; // index of the next vertex to make

    // Prepare number of statebits and priobits
    int statebits = 1;
    while ((1ULL<<(statebits)) <= (uint64_t)data->noStates) statebits++;
    int evenMax = 2 + 2*((data->noAccSets+1)/2); // should be enough...
    int priobits = 1;
    while ((1ULL<<(priobits)) <= (unsigned)evenMax) priobits++;

    std::vector<int> succ_state;  // for current state, the successors
    std::vector<int> succ_inter;  // for current intermediate state, the successors

    MTBDD trans_bdd = mtbdd_false;
    mtbdd_refs_pushptr(&trans_bdd);

    std::set<MTBDD> inter_bdds;
    std::set<uint64_t> targets;

    std::map<MTBDD, int> inter_vertices;
    std::map<uint64_t, int> target_vertices;

    // these variables are used deeper in the for loops, however we can
    // push them to mtbdd refs here and avoid unnecessary pushing and popping
    MTBDD lblbdd = mtbdd_false;
    MTBDD leaf = mtbdd_false;
    mtbdd_refs_pushptr(&lblbdd);
    mtbdd_refs_pushptr(&leaf);

    int ref_counter = 0;

    // Loop over every state
    for (int i=0; i<data->noStates; i++) {
        auto state = data->states+i;
        trans_bdd = mtbdd_false;

        // Loop over all transitions of the current state
        for (int j=0; j<state->noTrans; j++) {
            auto trans = state->transitions+j;
            // there should be a single successor per transition
            assert(trans->noSucc == 1);
            // there should be a label at state or transition level
            BTree* label;
            if (state->label != NULL) label = state->label;
            else label = trans->label;
            assert(label != NULL);
            // there should be a priority at state or transition level
            int priority = 0;
            if (state->accSig == NULL) {
                // there should be exactly one acceptance set!
                assert(trans->noAccSig == 1);
                // adjust priority
                priority = adjustPriority(trans->accSig[0], isMaxParity, controllerIsOdd, data->noAccSets);
            }
            // translate the label to a BDD
            lblbdd = CALL(evalLabel, label, data, variables);
            leaf = mtbdd_int64(((uint64_t)priority << 32) | (uint64_t)(trans->successors[0]));
            // add the transition to the transition BDD
            trans_bdd = mtbdd_ite(lblbdd, leaf, trans_bdd);
            lblbdd = leaf = mtbdd_false;
        }

        // At this point, we have the transitions from the state all in a neat single BDD.
        // Time to generate the split game fragment from the current state.

        collect_inter(trans_bdd, uap_count, inter_bdds);
        for (MTBDD inter_bdd : inter_bdds) {
            MTBDD targets_bdd = collect_targets(inter_bdd, targets, statebits, priobits);

            int vinter;
            auto it = inter_vertices.find(targets_bdd);
            if (it == inter_vertices.end()) {
                for (uint64_t lval : targets) {
                    int priority = (int)(lval >> 32);
                    int target = (int)(lval & 0xffffffff);

                    if (priority != 0) {
                        int vfin;
                        auto it = target_vertices.find(lval);
                        if (it == target_vertices.end()) {
                            vfin = nextIndex++;
                            game->init_vertex(vfin, priority, 0);
                            game->e_start(vfin);
                            game->e_add(vfin, target);
                            game->e_finish();
                            target_vertices.insert(std::make_pair(lval, vfin));
                        } else {
                            vfin = it->second;
                        }
                        succ_inter.push_back(vfin);
                    } else {
                        succ_inter.push_back(target);
                    }
                }

                vinter = nextIndex++;
                game->init_vertex(vinter, 0, 0, "from " + std::to_string(state->id));
                game->e_start(vinter);
                for (int to : succ_inter) game->e_add(vinter, to);
                game->e_finish();
                inter_vertices.insert(std::make_pair(targets_bdd, vinter));
                mtbdd_refs_push(targets_bdd);
                ref_counter++;
                succ_inter.clear();
            } else {
                vinter = it->second;
            }

            succ_state.push_back(vinter);
            targets.clear();
            target_vertices.clear();
        }

        // there should be a priority at state or transition level
        int priority;
        if (state->noAccSig != 0) {
            priority = adjustPriority(state->accSig[0], isMaxParity, controllerIsOdd, data->noAccSets);
        } else {
            priority = 0;
        }

        game->init_vertex(state->id, priority, 1, state->name ? state->name : std::to_string(state->id));
        game->e_start(state->id);
        for (int to : succ_state) game->e_add(state->id, to);
        game->e_finish();

        succ_state.clear();
        inter_bdds.clear();
        inter_vertices.clear();
    }

    // tell Oink we're done adding stuff, resize game to final size
    game->v_resize(nextIndex);

    mtbdd_refs_popptr(3);
    mtbdd_refs_pop(ref_counter);

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
            ("sym", "Solve the parity game using the internal symbolic solver")
            ("naive", "Use the naive splitting procedure (not recommended)")
            ("explicit", "Use the explicit splitting procedure (not recommended)")
            ("real", "Only check realizability (no synthesis)")
            ("bisim-game", "Apply bisimulation minimisation to the game")
            ("bisim-sol", "Apply bisimulation minimisation to the solution")
            ("bisim", "Apply bisimulation minimisation (--bisim-game and --bisim-sol)")
            ("onehot", "Use one-hot encoding for the states (recommended)")
            ("isop", "Convert BDDs to AIG using ISOP (instead of Shannon expansion)")
            ("compress", "Compress the generated AIG using ABC")
            ("drewrite", "Compress the generated AIG using ABCs commands drw and drf")
            ("best", "Try all combinations of bisim and isop and write the smallest AIG")
            ("no-solve", "Do not solve, halt after constructing the parity game")
            ("print-game", "Just print the parity game (implies no-solve)")
            ("print-witness", "Print the witness parity game")
            ("a,write-ascii", "Write ascii AIGER file")
            ("b,write-binary", "Write binary AIGER file")
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
    } catch (const cxxopts::exceptions::exception & e) {
        std::cout << "error parsing options: " << e.what() << std::endl;
        exit(0);
    }
}


VOID_TASK_0(gc_start)
{
    std::cerr << "starting garbage collection..." << std::endl;
}


VOID_TASK_0(gc_end)
{
    std::cerr << "garbage collection finished." << std::endl;
}


TASK_1(int, main_task, cxxopts::ParseResult*, _options)
{
    auto & options = *_options;

    bool verbose = options["verbose"].count() > 0;

    // First initialize the HOA data structure
    const double t_before_parsing = wctime();
    HoaData* data = (HoaData*)malloc(sizeof(HoaData));
    defaultsHoa(data);

    if (options.unmatched().size() == 0) {
        int ret = parseHoa(stdin, data);
        if (ret != 0) return ret;
    } else {
        std::string filename = options.unmatched()[0];
        FILE* f = fopen(filename.c_str(), "r");
        if (f == NULL) {
            std::cout << "file not found: " << filename << std::endl;
            return 0;
        }
        int ret = parseHoa(f, data);
        fclose(f);
        if (ret != 0) return ret;
    }

    const double t_after_parsing = wctime();
    if (verbose) {
        std::cerr << "\033[1;37mfinished parsing automaton in " << std::fixed << (t_after_parsing - t_before_parsing) << " sec.\033[m" << std::endl;
        std::cerr << "automaton has " << data->noStates << " states." << std::endl;
    }

    // First check if the automaton is a parity automaton
    bool isMaxParity = true;
    short controllerParity = 0;
    int ret = isParityGFG(data, &isMaxParity, &controllerParity);
    if (ret != 0) return ret;
    bool controllerIsOdd = controllerParity != 0;

    // Check if priorities are either all on state or all on transition, check state IDs for compatibility
    bool state_priorities = false;

    bool is_bad = false;
    for (int i=0; i<data->noStates; i++) {
        if (i != data->states[i].id) {
            std::cerr << "state " << i << " has an invalid id "  << data->states[i].id << "!" << std::endl;
            is_bad = true;
        }
        if (data->states[i].noAccSig != 0) {
            if (i == 0) {
                state_priorities = true;
            } else if (state_priorities == false) {
                std::cerr << "not every state has a priority!";
                return -1;
            }
        }
    }
    if (is_bad) return -1;

    if (verbose) {
        if (state_priorities) std::cerr << "priorities are on states" << std::endl;
        else std::cerr << "priorities are on transitions" << std::endl;
    }

    // Initialize Sylvan
    sylvan_set_limits(512LL << 22, 1, 14); // should be enough (2 gigabytes)
    sylvan_init_package();
    sylvan_init_mtbdd();
    sylvan_init_zdd();
    if (verbose) sylvan_gc_hook_pregc(TASK(gc_start));
    if (verbose) sylvan_gc_hook_postgc(TASK(gc_end));

    bool explicit_solver = options["sym"].count() == 0;
    bool naive_splitting = options["naive"].count() > 0;
    bool explicit_splitting = options["explicit"].count() > 0;
    bool write_pg = options["print-game"].count() > 0;
    bool no_solve = options["no-solve"].count() > 0;
    bool bisim_game = options["bisim"].count() > 0 or options["bisim-game"].count() > 0;
    bool bisim_sol = options["bisim"].count() > 0 or options["bisim-sol"].count() > 0;

    SymGame *sym = nullptr;
    bool realizable = false; // not known yet

    if (explicit_solver) {
        // Remember the start vertex
        int vstart = data->start[0];
        std::map<int, MTBDD> vertex_to_bdd;

        // Construct the explicit game
        pg::Game *game = nullptr;
        if (naive_splitting) {
            const double t_before_splitting = wctime();
            game = CALL(constructGameNaive, data, isMaxParity, controllerIsOdd);
            const double t_after_splitting = wctime();
            if (verbose) {
                std::cerr << "\033[1;37mfinished constructing game in " << std::fixed << (t_after_splitting - t_before_splitting) << " sec.\033[m" << std::endl;
            }
        } else if (explicit_splitting) {
            const double t_before_splitting = wctime();
            game = CALL(constructGame, data, isMaxParity, controllerIsOdd);
            const double t_after_splitting = wctime();
            if (verbose) {
                std::cerr << "\033[1;37mfinished constructing game in " << std::fixed << (t_after_splitting - t_before_splitting) << " sec.\033[m" << std::endl;
            }
        } else {
            const double t_1 = wctime();
            vstart = 0; // always set to 0 by constructSymGame
            sym = CALL(constructSymGame, data, isMaxParity, controllerIsOdd);
            const double t_2 = wctime();
            if (verbose) std::cerr << "\033[1;37mfinished constructing symbolic game in " << std::fixed << (t_2 - t_1) << " sec.\033[m" << std::endl;
            if (bisim_game) {
                const double t_before = wctime();
                MTBDD partition = CALL(min_lts_strong, sym, false);
                mtbdd_protect(&partition);
                CALL(minimize, sym, partition, verbose);
                mtbdd_unprotect(&partition);
                std::cerr << "number of blocks: " << count_blocks() << std::endl;
                const double t_after = wctime();
                if (verbose) std::cerr << "\033[1;37mfinished bisimulation minimisation of game in " << std::fixed << (t_after - t_before) << " sec.\033[m" << std::endl;
            }
            const double t_3 = wctime();
            game = sym->toExplicit(vertex_to_bdd);
            const double t_4 = wctime();
            if (verbose) {
                std::cerr << "\033[1;37mfinished constructing explicit game in " << std::fixed << (t_4 - t_3) << " sec.\033[m" << std::endl;
            }
        }

        if (verbose) {
            std::cerr << "constructed game has " << game->vertexcount() << " vertices and " << game->edgecount() << " edges." << std::endl;
        }

        if (write_pg) {
            // in case we want to write the file to PGsolver file format...
            game->set_label(vstart, "initial");
            game->write_pgsolver(std::cout);
            exit(0);
        }

        if (no_solve) {
            exit(0);
        }

        // we sort now, so we can track the initial state
        double begin = wctime();
        int *mapping = new int[game->vertexcount()];
        game->sort(mapping);

        // OK, fire up the engine
        std::stringstream log;

        std::string solver = "tl";
        pg::Solvers solvers;
        for (unsigned id=0; id<solvers.count(); id++) {
            if (options.count(solvers.label(id))) solver = solvers.label(id);
        }

        pg::Oink engine(*game, verbose ? std::cerr : log);
        engine.setTrace(0); //verbose ? 1 : 0); actually don´t -- maybe add a 2nd verbosity level later
        engine.setRenumber();
        engine.setSolver(solver);
        engine.setWorkers(-1);

        // and run the solver
        engine.run();
        double end = wctime();

        // report how long it all took
        if (verbose) {
            std::cerr << "\033[1;37mfinished solving game in " << std::fixed << (end - begin) << " sec.\033[m" << std::endl;
        }

        // undo the sorting ## FOR SOME REASON, THIS CAN BE SLOW!
        game->permute(mapping);
        delete[] mapping;

        realizable = game->winner[vstart] == 0;

        // finally, check if the initial vertex is won by controller or environment
        if (realizable) {
            if (sym != nullptr) {
                // now get the strategy from Oink...
                std::map<MTBDD, MTBDD> str;  // cap to priostate
                for (int v=0; v<game->vertexcount(); v++) {
                    // good cap states are owned and won by 0
                    if (game->owner(v) == 0 && game->winner[v] == 0) {
                        auto a = vertex_to_bdd[v];
                        auto b = vertex_to_bdd[game->strategy[v]];
                        str[a] = b;
                    }
                }

                if (!sym->applyStrategy(str)) {
                    std::cerr << "cannot apply strategy!" << std::endl;
                }
            }
        }
    } else {
        // Construct the game
        const double t_before_construct = wctime();
        sym = CALL(constructSymGame, data, isMaxParity, controllerIsOdd);
        const double t_after_construct = wctime();

        if (verbose) {
            std::cerr << "\033[1;37mfinished constructing symbolic game in " << std::fixed << (t_after_construct - t_before_construct) << " sec.\033[m" << std::endl;
        }

        if (bisim_game) {
            const double t_before = wctime();
            MTBDD partition = CALL(min_lts_strong, sym, false);
            mtbdd_protect(&partition);
            CALL(minimize, sym, partition, verbose);
            mtbdd_unprotect(&partition);
            std::cerr << "number of blocks: " << count_blocks() << std::endl;
            const double t_after = wctime();
            if (verbose) std::cerr << "\033[1;37mfinished bisimulation minimisation of game in " << std::fixed << (t_after - t_before) << " sec.\033[m" << std::endl;
        }

        if (write_pg) {
            // in case we want to write the file to PGsolver file format...
            std::map<int, MTBDD> vertex_to_bdd;
            auto pg = sym->toExplicit(vertex_to_bdd);
            pg->write_pgsolver(std::cout);
            exit(0);
        }

        if (no_solve) {
            exit(0);
        }

        const double t_before_solve = wctime();
        realizable = CALL(wrap_solve, sym, verbose);
        const double t_after_solve = wctime();

        if (verbose) {
            std::cerr << "\033[1;37mfinished solving game in " << std::fixed << (t_after_solve - t_before_solve) << " sec.\033[m" << std::endl;
        }
    }

    if (options["real"].count() > 0) {
        if (realizable) {
            std::cout << "REALIZABLE" << std::endl;
            if (verbose) {
                std::cerr << "\033[1;37mtotal time was " << std::fixed << (wctime() - t_before_parsing) << " sec.\033[m" << std::endl;
            }
            exit(10);
        } else {
            std::cout << "UNREALIZABLE" << std::endl;
            if (verbose) {
                std::cerr << "\033[1;37mtotal time was " << std::fixed << (wctime() - t_before_parsing) << " sec.\033[m" << std::endl;
            }
            exit(20);
        }
    }

    if (realizable) {
        if (verbose) std::cerr << "\033[1;38;5;10mgame is realizable!\033[m" << std::endl;

        if (naive_splitting or explicit_splitting) {
            std::cerr << "--naive and --explicit are incompatible with generating the AIG!" << std::endl;
            exit(10);
        }

        if (verbose) {
            // sym->print_trans();
            // sym->print_strategies();
        }

        const double t_before = wctime();
        CALL(wrap_pp, sym, verbose);
        const double t_after = wctime();

        if (verbose) {
            std::cerr << "\033[1;37mfinished post processing in " << std::fixed << (t_after - t_before) << " sec.\033[m" << std::endl;
        }

        if (verbose) {
            // sym->print_trans();
            // sym->print_strategies();
        }

        const bool best = options["best"].count() > 0;

        if (best) {
            AIGmaker var1(data, sym);
            var1.process();
            AIGmaker var2(data, sym);
            var2.setIsop();
            var2.process();
            AIGmaker var3(data, sym);
            var3.setOneHot();
            var3.process();

            const double t_before = wctime();
            MTBDD partition = RUN(min_lts_strong, sym, true);
            mtbdd_protect(&partition);
            RUN(minimize, sym, partition, verbose);
            mtbdd_unprotect(&partition);
            const double t_after = wctime();
            if (verbose) std::cerr << "\033[1;37mfinished bisimulation minimisation of solution in " << std::fixed << (t_after - t_before) << " sec.\033[m" << std::endl;

            AIGmaker var1b(data, sym);
            var1b.process();
            AIGmaker var2b(data, sym);
            var2b.setIsop();
            var2b.process();
            AIGmaker var3b(data, sym);
            var3b.setOneHot();
            var3b.process();

            if (verbose) {
                std::cerr << "no bisim, ite: " << var1.getNumAnds() << std::endl;
                std::cerr << "no bisim, isop: " << var2.getNumAnds() << std::endl;
                std::cerr << "no bisim, oh: " << var3.getNumAnds() << std::endl;
                std::cerr << "bisim, ite: " << var1b.getNumAnds() << std::endl;
                std::cerr << "bisim, isop: " << var2b.getNumAnds() << std::endl;
                std::cerr << "bisim, oh: " << var3b.getNumAnds() << std::endl;
            }

            if (options["drewrite"].count() > 0) {
                var1.drewrite();
                var2.drewrite();
                var3.drewrite();
                var1b.drewrite();
                var2b.drewrite();
                var3b.drewrite();

                if (verbose) {
                    std::cerr << "sizes after drw+drf with ABC:" << std::endl;
                    std::cerr << "no bisim, ite: " << var1.getNumAnds() << std::endl;
                    std::cerr << "no bisim, isop: " << var2.getNumAnds() << std::endl;
                    std::cerr << "no bisim, oh: " << var3.getNumAnds() << std::endl;
                    std::cerr << "bisim, ite: " << var1b.getNumAnds() << std::endl;
                    std::cerr << "bisim, isop: " << var2b.getNumAnds() << std::endl;
                    std::cerr << "bisim, oh: " << var3b.getNumAnds() << std::endl;
                }
            }

            if (options["compress"].count() > 0) {
                var1.compress();
                var2.compress();
                var3.compress();
                var1b.compress();
                var2b.compress();
                var3b.compress();

                if (verbose) {
                    std::cerr << "sizes after compressing with ABC:" << std::endl;
                    std::cerr << "no bisim, ite: " << var1.getNumAnds() << std::endl;
                    std::cerr << "no bisim, isop: " << var2.getNumAnds() << std::endl;
                    std::cerr << "no bisim, oh: " << var3.getNumAnds() << std::endl;
                    std::cerr << "bisim, ite: " << var1b.getNumAnds() << std::endl;
                    std::cerr << "bisim, isop: " << var2b.getNumAnds() << std::endl;
                    std::cerr << "bisim, oh: " << var3b.getNumAnds() << std::endl;
                }
            }

            auto smallest = var1.getNumAnds();
            smallest = std::min(smallest, var2.getNumAnds());
            smallest = std::min(smallest, var3.getNumAnds());
            smallest = std::min(smallest, var1b.getNumAnds());
            smallest = std::min(smallest, var2b.getNumAnds());
            smallest = std::min(smallest, var3b.getNumAnds());

            if (options.count("write-binary")) {
                if (var1.getNumAnds() == smallest) {
                    var1.writeBinary(stdout);
                } else if (var2.getNumAnds() == smallest) {
                    var2.writeBinary(stdout);
                } else if (var3.getNumAnds() == smallest) {
                    var3.writeBinary(stdout);
                } else if (var1b.getNumAnds() == smallest) {
                    var1b.writeBinary(stdout);
                } else if (var2b.getNumAnds() == smallest) {
                    var2b.writeBinary(stdout);
                } else if (var3b.getNumAnds() == smallest) {
                    var3b.writeBinary(stdout);
                }
            } else if (options["write-ascii"].count() > 0) {
                if (var1.getNumAnds() == smallest) {
                    var1.writeAscii(stdout);
                } else if (var2.getNumAnds() == smallest) {
                    var2.writeAscii(stdout);
                } else if (var3.getNumAnds() == smallest) {
                    var3.writeAscii(stdout);
                } else if (var1b.getNumAnds() == smallest) {
                    var1b.writeAscii(stdout);
                } else if (var2b.getNumAnds() == smallest) {
                    var2b.writeAscii(stdout);
                } else if (var3b.getNumAnds() == smallest) {
                    var3b.writeAscii(stdout);
                }
            }
            if (verbose) {
                std::cerr << "\033[1;37mtotal time was " << std::fixed << (wctime() - t_before_parsing) << " sec.\033[m" << std::endl;
            }
            exit(10);
        }

        if (bisim_sol) {
            const double t_before = wctime();
            MTBDD partition = CALL(min_lts_strong, sym, true);
            mtbdd_protect(&partition);
            if (verbose) {
                // CALL(print_partition, sym, partition);
            }
            CALL(minimize, sym, partition, verbose);
            mtbdd_unprotect(&partition);
            const double t_after = wctime();
            if (verbose) std::cerr << "\033[1;37mfinished bisimulation minimisation of solution in " << std::fixed << (t_after - t_before) << " sec.\033[m" << std::endl;
        }

        if (verbose) {
            // sym->print_trans();
            // sym->print_strategies();
        }

        /**
         * maybe print witness parity game, which should be fully won by even
         */

        if (options["print-witness"].count() > 0) {
            auto pargame = sym->strategyToPG();
            pargame->write_pgsolver(std::cout);
            if (verbose) {
                std::cerr << "\033[1;37mtotal time was " << std::fixed << (wctime() - t_before_parsing) << " sec.\033[m" << std::endl;
            }
            exit(10);
        }

        const double t_before_encoding = wctime();
        AIGmaker maker(data, sym);
        if (verbose) {
            maker.setVerbose();
        }
        if (options["isop"].count() > 0) {
            maker.setIsop();
        }
        if (options["onehot"].count() > 0) {
            maker.setOneHot();
        }
        {
            const double t_before = wctime();
            maker.process();
            const double t_after = wctime();
            if (verbose) std::cerr << "\033[1;37mfinished encoding in " << std::fixed << (t_after - t_before) << " sec.\033[m" << std::endl;
        }

        /**
         * maybe compress with ABC
         */
        if (options["drewrite"].count() > 0) {
            if (verbose) std::cerr << "size of AIG before drw+drf: " << maker.getNumAnds() << " gates." << std::endl;
            const double t_before = wctime();
            maker.drewrite();
            const double t_after = wctime();
            if (verbose) std::cerr << "size of AIG after drw+drf: " << maker.getNumAnds() << " gates." << std::endl;
            if (verbose) std::cerr << "\033[1;37mfinished drw+drf in " << std::fixed << (t_after - t_before) << " sec.\033[m" << std::endl;
        }

        if (options["compress"].count() > 0) {
            if (verbose) std::cerr << "size of AIG before compression: " << maker.getNumAnds() << " gates." << std::endl;
            const double t_before = wctime();
            maker.compress();
            const double t_after = wctime();
            if (verbose) std::cerr << "size of AIG after compression: " << maker.getNumAnds() << " gates." << std::endl;
            if (verbose) std::cerr << "\033[1;37mfinished compression in " << std::fixed << (t_after - t_before) << " sec.\033[m" << std::endl;
        }

        if (verbose) std::cerr << "final size of AIG: " << maker.getNumAnds() << " gates." << std::endl;


        if (options["write-binary"].count() > 0) {
            maker.writeBinary(stdout);
        } else if (options["write-ascii"].count() > 0) {
            maker.writeAscii(stdout);
        }
        if (verbose) {
            std::cerr << "\033[1;37mtotal time was " << std::fixed << (wctime() - t_before_parsing) << " sec.\033[m" << std::endl;
        }
        exit(10);
    } else {
        if (verbose) {
            std::cerr << "\033[1;37mtotal time was " << std::fixed << (wctime() - t_before_parsing) << " sec.\033[m" << std::endl;
            std::cerr << "\033[1;31mgame is unrealizable!\033[m" << std::endl;
        }
        exit(20);
    }

    if (verbose) {
        std::cerr << "\033[1;37mtotal time was " << std::fixed << (wctime() - t_before_parsing) << " sec.\033[m" << std::endl;
    }

    if (verbose) sylvan_stats_report(stdout);

    // We don't need Sylvan anymore at this point
    sylvan_quit();

    // free HOA allocated data structure
    resetHoa(data);
}



/**
 * The main function
 */
int
main(int argc, char* argv[])
{
    auto options = handleOptions(argc, argv);

    // Initialize Lace, only 1 worker
    lace_start(1, 1024*1024*2); // initialize Lace, but sequentially
                                // also get a large enough task size... (2M tasks) for PSI!

    int res = RUN(main_task, &options);

    lace_stop();

    return res;
}

