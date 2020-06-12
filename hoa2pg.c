/**************************************************************************
 * Copyright (c) 2019- Guillermo A. Perez
 * 
 * This file is part of HOATOOLS.
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
 *************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "simplehoa.h"

// A simple data structure for a linked list of PGSolver vertex
// entries

typedef struct PGSVertices {
    int id;
    int owner;
    int priority;
    char* name;
    IntList* successors;
    struct PGSVertices* next;
} PGSVertices;

static PGSVertices* initPGSVertices(size_t amount) {
    assert(amount > 0);
    PGSVertices* head = malloc(sizeof(PGSVertices));
    head->name = NULL;
    head->successors = NULL;
    head->next = NULL;
    PGSVertices* cur = head;
    for (int i = 0; i < amount - 1; i++) {
        cur->next = malloc(sizeof(PGSVertices));
        cur = cur->next;
        cur->name = NULL;
        cur->successors = NULL;
        cur->next = NULL;
    }
    return head;
}

static void deletePGSVertices(PGSVertices* u) {
    if (u == NULL)
        return;
    PGSVertices* v = u->next;
    if (u->name != NULL)
        free(u->name);
    IntList* succ = u->successors;
    while (succ != NULL) {
        IntList* temp = succ->next;
        free(succ);
        succ = temp;
    }
    free(u);
    deletePGSVertices(v);
}

static void printPGSVertices(PGSVertices* vlist) {
    for (PGSVertices* v = vlist; v != NULL; v = v->next) {
        assert(v->successors != NULL);
        printf("%d %d %d %d", v->id, v->priority, v->owner,
               v->successors->i);
        for (IntList* succ = v->successors->next; succ != NULL; succ = succ->next) {
            printf(",%d", succ->i);
        }
        if (v->name != NULL)
            printf(" \"%s\"\n", v->name);
        else
            printf(" \"%d\"\n", v->id);
    }
}

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
    fprintf(stderr, "Changed %d into %d. Original objective: %s %s with "
                    "maximal priority %d\n", p, shifted,
                    (maxPriority ? "max" : "min"),
                    (winRes == 0 ? "even" : "odd"),
                    noPriorities);
#endif
    return shifted;
}


/* Read the EHOA file, construct a graph-based game and
 * dump it as a PGSolver game
 */
int main(int argc, char* argv[]) {
    HoaData* data = malloc(sizeof(HoaData));
    defaultsHoa(data);
    int ret = parseHoa(stdin, data);
    // 0 means everything was parsed correctly
    if (ret != 0)
        return ret;
    // A few semantic checks! TODO: use checkParityGFG function instead
    // (1) the automaton should be a parity one
    if (strcmp(data->accNameID, "parity") != 0) {
        fprintf(stderr, "Expected \"parity...\" automaton, found \"%s\" "
                        "as automaton type\n", data->accNameID);
        return 100;
    }
    bool foundOrd = false;
    bool maxPriority;
    bool foundRes = false;
    short winRes;
    for (StringList* param = data->accNameParameters; param != NULL;
            param = param->next) {
        if (strcmp(param->str, "max") == 0) {
            maxPriority = true;
            foundOrd = true;
        }
        if (strcmp(param->str, "min") == 0) {
            maxPriority = false;
            foundOrd = true;
        }
        if (strcmp(param->str, "even") == 0) {
            winRes = 0;
            foundRes = true;
        }
        if (strcmp(param->str, "odd") == 0) {
            winRes = 1;
            foundRes = true;
        }
    }
    if (!foundOrd) {
        fprintf(stderr, "Expected \"max\" or \"min\" in the acceptance name\n");
        return 101;
    }
    if (!foundRes) {
        fprintf(stderr, "Expected \"even\" or \"odd\" in the acceptance name\n");
        return 102;
    }
    // (2) the automaton should be deterministic, complete, colored
    bool det = false;
    bool complete = false;
    bool colored = false;
    for (StringList* prop = data->properties; prop != NULL; prop = prop->next) {
        if (strcmp(prop->str, "deterministic") == 0)
            det = true;
        if (strcmp(prop->str, "complete") == 0)
            complete = true;
        if (strcmp(prop->str, "colored") == 0)
            colored = true;
    }
    if (!det) {
        fprintf(stderr, "Expected a deterministic automaton, "
                        "did not find \"deterministic\" in the properties\n");
        return 200;
    }
    if (!complete) {
        fprintf(stderr, "Expected a complete automaton, "
                        "did not find \"complete\" in the properties\n");
        return 201;
    }
    if (!colored) {
        fprintf(stderr, "Expected one acceptance set per transition, "
                        "did not find \"colored\" in the properties\n");
        return 202;
    }
    // (3) the automaton should have a unique start state
    if (data->start == NULL || data->start->next != NULL) {
        fprintf(stderr, "Expected a unique start state\n");
        return 300;
    }

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
            fprintf(stderr, "Found an uncontrollable AP: %d\n", i);
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
    int nextIndex = data->noStates;
    // we need to store vertices and edges somehow:
    PGSVertices* states = initPGSVertices(data->noStates);
    PGSVertices* curState = states;
    PGSVertices* partVals = initPGSVertices(data->noStates * numValuations);
    PGSVertices* curPartVal = partVals;
    PGSVertices* fullVals = NULL;
    PGSVertices* curFullVal = NULL;
    for (StateList* state = data->states; state != NULL;
         state = state->next) {
        int firstSucc = nextIndex;
        nextIndex += numValuations;
        for (unsigned value = 0; value < numValuations; value++) {
            IntList* validVals = NULL;
            for (TransList* trans = state->transitions;
                 trans != NULL; trans = trans->next) {
                // there should be a single successor per transition
                assert(trans->successors != NULL &&
                       trans->successors->next == NULL);
                // there should be a label at state or transition level
                BTree* label;
                if (state->label != NULL)
                    label = state->label;
                else
                    label = trans->label;
                assert(label != NULL);
                // there should be a priority at state or transition level
                IntList* acc = state->accSig;
                if (state->accSig == NULL)
                    acc = trans->accSig;
                // one of the two should be non-NULL
                // and there should be exactly one acceptance set!
                assert(acc != NULL && acc->next == NULL);
                int priority = adjustPriority(acc->i, maxPriority, winRes,
                                              data->noAccSets);
                // we add a vertex + edges if the transition is compatible with the
                // valuation we are currently considering; because of PGSolver
                // format we add only the leaving edge to a state and defer
                // edges to it (from partial valuations) to later
                int evald = evalLabel(label, data->aliases, numUcntAPs,
                                      ucntAPs, value);
#ifndef NDEBUG
                fprintf(stderr, "Called evalLabel for value %d with "
                                "%d uncontrollable APs; got %d\n",
                        value, numUcntAPs, evald);
#endif
                if (evald != -1) {
                    int fval = nextIndex;
                    nextIndex++;
                    // as unique successor we add the successor via the
                    // transition (so the choice of player is unimportant)
                    if (curFullVal == NULL) {
                        assert(fullVals == NULL);
                        fullVals = malloc(sizeof(PGSVertices));
                        curFullVal = fullVals;
                    } else {
                        assert(fullVals != NULL);
                        curFullVal->next = malloc(sizeof(PGSVertices));
                        curFullVal = curFullVal->next;
                    }
                    assert(curFullVal != NULL);
                    curFullVal->id = fval;
                    curFullVal->priority = priority;
                    curFullVal->owner = 0;
                    curFullVal->name = NULL;
                    curFullVal->next = NULL;
                    curFullVal->successors = prependIntNode(NULL, trans->successors->i);
                    // we also keep it as a list for the successors of the
                    // partial value
                    validVals = prependIntNode(validVals, fval);
                }
            }
            assert(validVals != NULL);
            assert(curPartVal != NULL);
            // Now we can add priority-0 edges from the player-0
            // partial-valuation vertices to full-valuation vertices
            curPartVal->id = firstSucc + value;
            curPartVal->priority = 0;
            curPartVal->owner = 0;
            curPartVal->successors = validVals;
            curPartVal = curPartVal->next;
        }
        // Now we can add priority-0 edges from the player-1 vertex to
        // all partial-valuation vertices owned by player 0
        assert(curState != NULL);
        curState->id = state->id;
        curState->priority = 0;
        curState->owner = 1;
        for (int succ = firstSucc; succ < firstSucc + numValuations; succ++) 
            curState->successors = prependIntNode(curState->successors, succ);
        curState->name = state->name;
        curState = curState->next;
    }

    // Step 3: print the PGSolver file
    printf("parity %d;\n", nextIndex - 1);
    printPGSVertices(states);
#ifndef NDEBUG
    fprintf(stderr, "Partial valuation vertices:\n");
#endif
    printPGSVertices(partVals);
#ifndef NDEBUG
    fprintf(stderr, "Full valuation vertices:\n");
#endif
    printPGSVertices(fullVals);

    // Free dynamic memory
    deleteHoa(data);
    deletePGSVertices(states);
    deletePGSVertices(partVals);
    deletePGSVertices(fullVals);
    return EXIT_SUCCESS;
}
