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

extern "C" {
#include "simplehoa.h"
}

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
static int evalLabel(BTree* label, AliasList* aliases, 
                      int numAPs, int* apIds, unsigned value) {
    assert(label != NULL);
    int left;
    int right;
    unsigned mask;
    switch (label->type) {
        case NT_BOOL:
            return label->id ? 1 : -1;  // 0 becomes -1 like this
        case NT_AND:
            left = evalLabel(label->left, aliases, numAPs, apIds, value);
            right = evalLabel(label->right, aliases, numAPs, apIds, value);
            if (left == -1 || right == -1)
                return -1;
            if (left == 0 || right == 0)
                return 0;
            // otherwise
            return 1;
        case NT_OR:
            left = evalLabel(label->left, aliases, numAPs, apIds, value);
            right = evalLabel(label->right, aliases, numAPs, apIds, value);
            if (left == 1 || right == 1)
                return 1;
            if (left == 0 || right == 0)
                return 0;
            // otherwise
            return -1;
        case NT_NOT:
            return -1 * evalLabel(label->left, aliases, numAPs, apIds, value);
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
            for (AliasList* a = aliases; a != NULL; a = a->next) {
                if (strcmp(a->alias, label->alias) == 0)
                    return evalLabel(a->labelExpr, aliases, numAPs, apIds, value);
            }
            break;
        default:
            assert(false);  // all cases should be covered above
    }
    return -2;
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

    // Step 1: prepare a list of all uncontrollable inputs
    int numUcntAPs = data->noAPs;
    for (IntList* c = data->cntAPs; c != NULL; c = c->next)
        numUcntAPs--;
    int ucntAPs[numUcntAPs];
    int apIdx = 0;
    for (int i = 0; i < data->noAPs; i++) {
        bool found = false;
        for (IntList* c = data->cntAPs; c != NULL; c = c->next) {
            if (i == c->i) {
                found = true;
                break;
            }
        }
        if (!found) {
            ucntAPs[apIdx] = i;
#ifndef NDEBUG
            if (0) {
                fprintf(stderr, "Found an uncontrollable AP: %d\n", i);
            }
#endif
            apIdx++;
        }
    }
    assert(apIdx == numUcntAPs);

    // Step 2: for all states in the automaton and all valuations, create
    // vertices for both players and edges to go with them
    // NOTE: states retain their index while "intermediate" state-valuation
    // vertices receive new indices
    const unsigned numValuations = (1 << numUcntAPs);
    int nextIndex = data->noStates; // number of states

    int nv = data->noStates + data->noStates*numValuations; // it will actually grow larger...
    pg::Game game(nv);

    std::vector<int> inter;
    std::vector<int> succie;

    // Loop over every state
    for (StateList* state = data->states; state != NULL; state = state->next) {
        for (unsigned value = 0; value < numValuations; value++) {
            // for every valuation to the uncontrollable APs, we make an intermediate vertex
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
                // we add a vertex + edges if the transition is compatible with the
                // valuation we are currently considering
                int evald = evalLabel(label, data->aliases, numUcntAPs, ucntAPs, value);
                if (evald == -1) continue; // not compatible
                // adjust priority
                int priority = adjustPriority(acc->i, maxPriority, winRes, data->noAccSets);

                int vfin = nextIndex++;
                game.init_vertex(vfin, priority, 0);
                game.e_start(vfin);
                game.e_add(vfin, trans->successors->i);
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
        }

        game.init_vertex(state->id, 0, 1, state->name ? state->name : "");
        game.e_start(state->id);
        for (int to : inter) game.e_add(state->id, to);
        game.e_finish();
        inter.clear();
    }

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
        return 10;
    } else {
        std::cout << "UNREALIZABLE";
        return 20;
    }
}
