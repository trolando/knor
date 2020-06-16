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

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sys/time.h>

#include <game.hpp>
#include <oink.hpp>
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

static double t_start;

/* Given a label and a valuation of some of the atomic propositions,
 * we determine whether the label is true (1), false (-1), or its
 * value is unknown (0). The valuation is expected as an unsigned
 * integer whose i-th bit is 1 iff the i-th AP in apIds is set to 1
 */
#define evalLabel(a,b,c) CALL(evalLabel, a, b, c)
TASK_3(MTBDD, evalLabel, BTree*, label, AliasList*, aliases, uint32_t*, variables)
{
    assert(label != NULL);
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

/* Adjust the priorities since we have to make sure we output a max, even
 * parity game and that the priorities of player-0 vertices are useless
 * (that is, irrelevant). This assumes maxPriority is true iff the original
 * objective is max and winRes = 0 if even, otherwise it is 1 if odd.
 */
static inline int adjustPriority(int p, bool maxPriority, short winRes,
                                 int noPriorities) {
    // To deal with max vs min, we subtract from noPriorities if
    // originally it was min (for this we need it to be even!)
    int evenMax = noPriorities;
    if (evenMax % 2 != 0)
        evenMax += 1;
    int pForMax = maxPriority ? p : evenMax - p;
    // The plan is to use 0 as the priority for player-0 vertices,
    // this means shifting everything up; we take the opportunity
    // to make odd priorities even if the original objective asked
    // for odd ones
    int shifted = pForMax + (2 - winRes);
#ifndef NDEBUG
    if (0) {
        fprintf(stderr, "Changed %d into %d. Original objective: %s %s with "
                    "maximal priority %d\n", p, shifted,
                    (maxPriority ? "max" : "min"),
                    (winRes == 0 ? "even" : "odd"),
                    noPriorities);
    }
#endif
    return shifted;
}



void
collect_inter(MTBDD trans, int uap_count, int cap_count, std::set<MTBDD> &res)
{
    if (mtbdd_isleaf(trans)) {
        res.insert(trans);
    } else {
        // node
        uint32_t var = mtbdd_getvar(trans);
        if (var < (unsigned)uap_count) {
            // uncontrollable ap
            collect_inter(mtbdd_gethigh(trans), uap_count, cap_count, res);
            collect_inter(mtbdd_getlow(trans), uap_count, cap_count, res);
        } else {
            // controllable ap
            res.insert(trans);
        }
    }
}


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



int
main(int argc, char* argv[])
{
    t_start = wctime();

    // First initialize the HOA data structure
    HoaData* data = (HoaData*)malloc(sizeof(HoaData));
    defaultsHoa(data);

    if (argc == 1) {
        int ret = parseHoa(stdin, data);
        if (ret != 0) return ret;
    } else {
        FILE* f = fopen(argv[1], "r");
        int ret = parseHoa(f, data);
        fclose(f);
        if (ret != 0) return ret;
    }

    std::cerr << "finished reading file." << std::endl;

    // First check if the automaton is a parity automaton
    bool maxPriority = true;
    short int winRes = 0;

    int ret = isParityGFG(data, &maxPriority, &winRes);
    if (ret != 0) return ret;

    // Step 1: which APs are controllable
    pg::bitset controllable(data->noAPs);
    for (IntList* c = data->cntAPs; c != NULL; c = c->next) controllable[c->i] = 1;

    const int cap_count = controllable.count();
    const int uap_count = data->noAPs - cap_count;

    // order uncontrollable < controllable
    uint32_t variables[data->noAPs];
    long uidx = 0, oidx = data->noAPs - controllable.count();
    for (int i=0; i<data->noAPs; i++) {
        if (controllable[i]) variables[i] = oidx++;
        else variables[i] = uidx++;
    }
    assert(uidx == (data->noAPs - controllable.count()));
    assert(oidx == data->noAPs);

    // At this point, we need to initialize Sylvan for support
    lace_init(1, 0); // initialize Lace, but sequentially
    lace_startup(0, 0, 0); // no thread spawning
    LACE_ME;
    sylvan_set_limits(16LL << 25, 1, 16); // should be enough
    sylvan_init_package();
    sylvan_init_mtbdd();

    int nv = data->noStates * 100; // it will actually grow larger automatically
    pg::Game game(nv);
    int nextIndex = data->noStates; // number of states

    std::vector<int> inter;
    std::vector<int> succie;

    MTBDD var_bdd = mtbdd_set_from_array(variables, data->noAPs);
    mtbdd_refs_pushptr(&var_bdd);

    MTBDD uvar_bdd = mtbdd_set_from_array(variables, uap_count);
    mtbdd_refs_pushptr(&uvar_bdd);

    MTBDD cvar_bdd = mtbdd_set_from_array(variables+uap_count, cap_count);
    mtbdd_refs_pushptr(&cvar_bdd);

    MTBDD trans_bdd = mtbdd_false;
    mtbdd_refs_pushptr(&trans_bdd);

    std::set<MTBDD> halfway_bdds;
    std::set<uint64_t> targets;

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
            IntList* acc = state->accSig;
            if (state->accSig == NULL) acc = trans->accSig;
            // there should be exactly one acceptance set!
            assert(acc != NULL && acc->next == NULL);

            // adjust priority
            int priority = adjustPriority(acc->i, maxPriority, winRes, data->noAccSets);
            int target = trans->successors->i;

            // tricky tricky
            MTBDD lblbdd = evalLabel(label, data->aliases, variables);
            MTBDD leaf = mtbdd_int64(((uint64_t)priority << 32) | (uint64_t)target);
            trans_bdd = mtbdd_ite(lblbdd, leaf, trans_bdd);
        }

        collect_inter(trans_bdd, uap_count, cap_count, halfway_bdds);

        for (MTBDD ASDF : halfway_bdds) {
            collect_targets(ASDF, targets);
            for (uint64_t lval : targets) {
                int priority = (int)(lval >> 32);
                int target = (int)(lval & 0xffffffff);

                int vfin = nextIndex++;
                game.init_vertex(vfin, priority, 0);
                game.e_start(vfin);
                game.e_add(vfin, target);
                game.e_finish();
                succie.push_back(vfin);
            }

            int vinter = nextIndex++;
            inter.push_back(vinter);
            game.init_vertex(vinter, 0, 0);
            game.e_start(vinter);
            for (int to : succie) game.e_add(vinter, to);
            game.e_finish();
            succie.clear();
            targets.clear();
        }

        game.init_vertex(state->id, 0, 1, state->name ? state->name : "");
        game.e_start(state->id);
        for (int to : inter) game.e_add(state->id, to);
        game.e_finish();
        inter.clear();
        halfway_bdds.clear();
    }

    mtbdd_refs_popptr(4);
    sylvan_quit();

    int vstart = data->start->i;

    deleteHoa(data);
    std::cerr << "finished constructing game." << std::endl;

    game.v_resize(nextIndex);

    // game.write_pgsolver(std::cout);
    // std::cerr << "initial vertex: " << vstart << std::endl;

    int *mapping = new int[game.vertexcount()];
    game.sort(mapping);
    for (int i=0; i<nextIndex; i++) {
        if (mapping[i] == vstart) {
            vstart = i;
            break;
        }
    }
    delete[] mapping;

    pg::Oink en(game, std::cerr);
    en.setTrace(0);
    en.setRenumber();
    if (argc > 2) en.setSolver(argv[2]);
    else en.setSolver("fpi");
    en.setWorkers(-1);

    double begin = wctime();
    en.run();
    double end = wctime();

    std::cerr << "total solving time: " << std::fixed << (end-begin) << " sec." << std::endl;

    // game.write_pgsolver(std::cout);
    // std::cerr << "initial vertex: " << vstart << std::endl;
    // game.write_sol(std::cout);

    if (game.winner[vstart] == 0) {
        std::cout << "REALIZABLE";
        exit(10);
    } else {
        std::cout << "UNREALIZABLE";
        exit(20);
    }
}
