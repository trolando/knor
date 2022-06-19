/**************************************************************************
 * Copyright (c) 2019- Guillermo A. Perez
 * Copyright (c) 2020-2021 Tom van Dijk
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
#include <map>
#include <set>
#include <sys/time.h> // for gettimeofday

#include <game.hpp>
#include <oink.hpp>
#include <solvers.hpp>
#include <tools/cxxopts.hpp>
#include <sylvan.h>

extern "C" {
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
 * evalLabel converts a label on a transition to a MTBDD encoding the label
 * a label is essentially a boolean combination of atomic propositions, and aliases (aliases are rare)
 */
#define evalLabel(a,b,c) RUN(evalLabel,a,b,c)
TASK_3(MTBDD, evalLabel, BTree*, label, HoaData*, data, uint32_t*, variables)
{
    MTBDD left = mtbdd_false, right = mtbdd_false, result = mtbdd_false;
    mtbdd_refs_pushptr(&left);
    mtbdd_refs_pushptr(&right);
    switch (label->type) {
        case NT_BOOL:
            result = label->id ? mtbdd_true : mtbdd_false;
            break;
        case NT_AND:
            left = CALL(evalLabel, label->left, data, variables);
            right = CALL(evalLabel, label->right, data, variables);
            result = sylvan_and(left, right);
            break;
        case NT_OR:
            left = CALL(evalLabel, label->left, data, variables);
            right = CALL(evalLabel, label->right, data, variables);
            result = sylvan_or(left, right);
            break;
        case NT_NOT:
            left = CALL(evalLabel, label->left, data, variables);
            result = sylvan_not(left);
            break;
        case NT_AP:
            result = mtbdd_ithvar(variables[label->id]);
            break;
        case NT_ALIAS:
            // apply the alias
            for (int i=0; i<data->noAliases; i++) {
                Alias *a = data->aliases+i;
                if (strcmp(a->alias, label->alias) == 0) {
                    return CALL(evalLabel, a->labelExpr, data, variables);
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
 * Given a label and a valuation of some of the atomic propositions,
 * we determine whether the label is true (1), false (-1), or its
 * value is unknown (0). The valuation is expected as an unsigned
 * integer whose i-th bit is 1 iff the i-th AP in apIds is set to 1
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
 * Given some intermediary MTBDD root, collect all the MTBDD
 * leafs, representing the target vertices of the full transition
 * and the transition priority encoded in a single 64-bit value
 *
 * The input is essentially the result of collect_inter, i.e., just the controllable APs
 */
/*void
collect_targets(MTBDD trans, std::set<uint64_t> &res)
{
    if (mtbdd_isleaf(trans)) {
        res.insert(mtbdd_getint64(trans));
    } else {
        collect_targets(mtbdd_gethigh(trans), res);
        collect_targets(mtbdd_getlow(trans), res);
    }
}*/



MTBDD
encode_state(uint32_t state, const int statebits, const int priobits, const int offset)
{
    // create a cube
    MTBDD cube = mtbdd_true;
    for (int i=0; i<statebits; i++) {
        if (state & 1) cube = mtbdd_makenode(offset+statebits+priobits-i-1, mtbdd_false, cube);
        else cube = mtbdd_makenode(offset+statebits+priobits-i-1, cube, mtbdd_false);
        state >>= 1;
    }
    return cube;
}



MTBDD
encode_priostate(uint32_t state, uint32_t priority, const int statebits, const int priobits, const int offset)
{
    // create a cube
    MTBDD cube = mtbdd_true;
    for (int i=0; i<statebits; i++) {
        if (state & 1) cube = mtbdd_makenode(offset+statebits+priobits-i-1, mtbdd_false, cube);
        else cube = mtbdd_makenode(offset+statebits+priobits-i-1, cube, mtbdd_false);
        state >>= 1;
    }
    for (int i=0; i<priobits; i++) {
        if (priority & 1) cube = mtbdd_makenode(offset+priobits-i-1, mtbdd_false, cube);
        else cube = mtbdd_makenode(offset+priobits-i-1, cube, mtbdd_false);
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
collect_targets2(MTBDD trans, std::set<uint64_t> &res, const int statebits, const int priobits)
{
    if (mtbdd_isleaf(trans)) {
        uint64_t leaf = mtbdd_getint64(trans);
        res.insert(leaf);

        uint32_t priority = (uint32_t)(leaf>>32);
        uint32_t state = (uint32_t)(leaf & 0xffffffff);

        // std::cerr << "encoding " << priority << " " << state << std::endl;

        // create a cube
        MTBDD cube = mtbdd_true;
        for (int i=0; i<statebits; i++) {
            if (state & 1) cube = mtbdd_makenode(statebits+priobits-i-1, mtbdd_false, cube);
            else cube = mtbdd_makenode(statebits+priobits-i-1, cube, mtbdd_false);
            state >>= 1;
        }
        for (int i=0; i<priobits; i++) {
            if (priority & 1) cube = mtbdd_makenode(priobits-i-1, mtbdd_false, cube);
            else cube = mtbdd_makenode(priobits-i-1, cube, mtbdd_false);
            priority >>= 1;
        }

        return cube;
    } else {
        MTBDD left = collect_targets2(mtbdd_getlow(trans), res, statebits, priobits);
        mtbdd_refs_push(left);
        MTBDD right = collect_targets2(mtbdd_gethigh(trans), res, statebits, priobits);
        mtbdd_refs_push(right);
        MTBDD res = sylvan_or(left, right);
        mtbdd_refs_pop(2);
        return res;
    }
}



/**
 * Construct and solve the game explicitly
 */
pg::Game*
constructGameNaive(HoaData *data, bool isMaxParity, bool controllerIsOdd)
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
pg::Game*
constructGame(HoaData *data, bool isMaxParity, bool controllerIsOdd)
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
    // notice that the number of vertices automatically grows when needed anyway

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
            assert(trans->noSucc == 1);
            // tricky tricky
            lblbdd = evalLabel(label, data, variables);
            leaf = mtbdd_int64(((uint64_t)priority << 32) | (uint64_t)(trans->successors[0]));
            trans_bdd = mtbdd_ite(lblbdd, leaf, trans_bdd);
            lblbdd = leaf = mtbdd_false;
        }

        // At this point, we have the transitions from the state all in a neat
        // single BDD. Time to generate the split game fragment from the current
        // state.

        collect_inter(trans_bdd, uap_count, inter_bdds);
        for (MTBDD inter_bdd : inter_bdds) {
            MTBDD targets_bdd = collect_targets2(inter_bdd, targets, statebits, priobits);

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
                game->init_vertex(vinter, 0, 0);
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
    }

    // tell Oink we're done adding stuff, resize game to final size
    game->v_resize(nextIndex);

    mtbdd_refs_popptr(3);
    mtbdd_refs_pop(ref_counter);

    return game;
}


class SymGame {
public:
    int maxprio;
    int statebits;
    int priobits;
    int cap_count;
    int uap_count;

    MTBDD trans;       // transition relation of symbolic game:        state -> uap -> cap -> priority -> next_state
    MTBDD strategies;  // contains the solution: the strategies:       good_state -> uap -> cap
    MTBDD* cap_bdds;   // contains the solution: controllable ap bdds: state -> uap -> B
    MTBDD* state_bdds; // contains the solution: state bit bdds      : state -> uap -> B

    SymGame(int statebits, int priobits, int uap_count, int cap_count, int maxprio) {
        this->maxprio = maxprio;
        this->uap_count = uap_count;
        this->cap_count = cap_count;
        this->statebits = statebits;
        this->priobits = priobits;

        trans = mtbdd_false;
        strategies = mtbdd_false;
        mtbdd_protect(&trans);
        mtbdd_protect(&strategies);

        cap_bdds = new MTBDD[cap_count];
        for (int i=0; i<cap_count; i++) cap_bdds[i] = mtbdd_false;
        for (int i=0; i<cap_count; i++) mtbdd_protect(&cap_bdds[i]);
        state_bdds = new MTBDD[statebits];
        for (int i=0; i<statebits; i++) state_bdds[i] = mtbdd_false;
        for (int i=0; i<statebits; i++) mtbdd_protect(&state_bdds[i]);
    }

    virtual ~SymGame() {
        for (int i=0; i<cap_count; i++) mtbdd_unprotect(&cap_bdds[i]);
        for (int i=0; i<statebits; i++) mtbdd_unprotect(&state_bdds[i]);

        delete[] cap_bdds;
        delete[] state_bdds;
    }
};



/**
 * Construct the symbolic game
 */
SymGame*
constructSymGame(HoaData *data, bool isMaxParity, bool controllerIsOdd)
{
    int vstart = data->start[0];

    // Set which APs are controllable in the bitset controllable
    pg::bitset controllable(data->noAPs);
    for (int i=0; i<data->noCntAPs; i++) {
        controllable[data->cntAPs[i]] = 1;
    }

    // Count the number of controllable/uncontrollable APs
    const int cap_count = controllable.count();
    const int uap_count = data->noAPs - cap_count;

    // Prepare number of statebits and priobits
    int statebits = 1;
    while ((1ULL<<(statebits)) <= (uint64_t)data->noStates) statebits++;

    const int evenMax = 2 + 2*((data->noAccSets+1)/2); // should be enough...
    int priobits = 1;
    while ((1ULL<<(priobits)) <= (unsigned)evenMax) priobits++;

    SymGame *res = new SymGame(statebits, priobits, uap_count, cap_count, 0);

    // first <priobits, statebits> will be state variables

    // Initialize the BDD variable indices
    // Variable order uncontrollable < controllable
    uint32_t variables[data->noAPs];
    uint32_t uidx = statebits + priobits;
    uint32_t oidx = statebits + priobits + data->noAPs - controllable.count();
    for (int i=0; i<data->noAPs; i++) {
        if (controllable[i]) variables[i] = oidx++;
        else variables[i] = uidx++;
    }

    // only if verbose: std::cout << "prio bits: " << priobits << " and state bits: " << statebits << std::endl;

    // these variables are used deeper in the for loops, however we can
    // push them to mtbdd refs here and avoid unnecessary pushing and popping
    MTBDD transbdd = mtbdd_false;
    MTBDD statebdd = mtbdd_false;
    MTBDD lblbdd = mtbdd_false;
    MTBDD leaf = mtbdd_false;

    mtbdd_refs_pushptr(&transbdd);
    mtbdd_refs_pushptr(&statebdd);
    mtbdd_refs_pushptr(&lblbdd);
    mtbdd_refs_pushptr(&leaf);

    // Loop over every state
    for (int i=0; i<data->noStates; i++) {
        auto state = data->states+i;

        for (int j=0; j<state->noTrans; j++) {
            auto trans = state->transitions+j;

            // there should be a single successor per transition
            assert(trans->noSucc == 1);
            // there should be a label at state or transition level
            BTree* label = state->label != NULL ? state->label : trans->label;
            assert(label != NULL);
            // there should be a priority at state or transition level
            int priority = 0;
            if (trans->noAccSig != 0) {
                // there should be exactly one acceptance set!
                assert(trans->noAccSig == 1);
                // adjust priority
                priority = adjustPriority(trans->accSig[0], isMaxParity, controllerIsOdd, data->noAccSets);
            } else {
                auto id = trans->successors[0];
                assert(data->states[id].noAccSig == 1);
                priority = adjustPriority(data->states[id].accSig[0], isMaxParity, controllerIsOdd, data->noAccSets);
            }
            if (priority > res->maxprio) res->maxprio = priority;
            // swap with initial if needed
            int succ = trans->successors[0];
            if (succ == 0) succ = vstart;
            else if (succ == vstart) succ = 0;
            // encode the label as a MTBDD
            lblbdd = evalLabel(label, data, variables);
            // encode priostate (leaf) and update transition relation
            leaf = encode_priostate(succ, priority, statebits, priobits, statebits+priobits+data->noAPs);
            // trans := lbl THEN leaf ELSE trans
            transbdd = mtbdd_ite(lblbdd, leaf, transbdd);
            // deref lblbdd and leaf
            lblbdd = leaf = mtbdd_false;

            // std::cout << "added transition from " << state->id << " with prio " << priority << " to " << trans->successors[0] << std::endl; 
        }

        // encode source state and add to full transition relation
        int src = state->id;
        if (src == 0) src = vstart;
        else if (src == vstart) src = 0;
        statebdd = encode_state(src, statebits, priobits, 0);
        // update full trans with <statebdd> then <transbdd>
        res->trans = mtbdd_ite(statebdd, transbdd, res->trans);
        // deref statebdd and transbdd
        statebdd = transbdd = mtbdd_false;
    }

    mtbdd_refs_popptr(4);

    return res;
}


/**
 * Encode priority i.e. all states via priority <priority>
 */
MTBDD
encodeprio(int priority, int priobits)
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
 * Given a strategy, if there is choice, choose the zero edge
 */
TASK_2(MTBDD, clarify, MTBDD, str, MTBDD, cap_vars)
{
    if (mtbdd_set_isempty(cap_vars)) return str;

    uint32_t next_cap = mtbdd_getvar(cap_vars);
    if (mtbdd_isleaf(str)) {
        // str is leaf but we have variables to set to 0!
        MTBDD low = CALL(clarify, str, mtbdd_set_next(cap_vars));
        return mtbdd_makenode(next_cap, low, mtbdd_false);
    } else {
        uint32_t var = mtbdd_getvar(str);
        if (var < next_cap) {
            // not at controllable variables yet!
            MTBDD low = CALL(clarify, mtbdd_getlow(str), cap_vars);
            mtbdd_refs_pushptr(&low);
            MTBDD high = CALL(clarify, mtbdd_gethigh(str), cap_vars);
            mtbdd_refs_popptr(1);
            return mtbdd_makenode(var, low, high);
        } else if (next_cap < var) {
            // at controllable variables, and skipping one! (set to 0)
            MTBDD low = CALL(clarify, str, mtbdd_set_next(cap_vars));
            return mtbdd_makenode(next_cap, low, mtbdd_false);
        } else {
            // at controllable variables, and matching!
            if (mtbdd_getlow(str) != mtbdd_false) {
                // we can take the low branch, so only take the low branch
                MTBDD low = CALL(clarify, mtbdd_getlow(str), mtbdd_set_next(cap_vars));
                return mtbdd_makenode(next_cap, low, mtbdd_false);
            } else {
                // unfortunately, we cannot take the low branch!
                MTBDD high = CALL(clarify, mtbdd_gethigh(str), mtbdd_set_next(cap_vars));
                return mtbdd_makenode(next_cap, mtbdd_false, high);
            }
        }
    }
}


pg::Game*
symgame_to_pg(const SymGame &game)
{
    // Construct s_vars the set of state variables of the symbolic game
    MTBDD s_vars = mtbdd_set_empty();
    mtbdd_refs_pushptr(&s_vars);
    for (int i=0; i<game.statebits; i++) s_vars = mtbdd_set_add(s_vars, game.priobits+i);

    // Construct uap_vars the set of uncontrollable AP variables
    MTBDD uap_vars = mtbdd_set_empty();
    mtbdd_refs_pushptr(&uap_vars);
    for (int i=0; i<game.uap_count; i++) uap_vars = mtbdd_set_add(uap_vars, game.priobits+game.statebits+i);

    // Construct cap_vars the set of controllable AP variables
    MTBDD cap_vars = mtbdd_set_empty();
    mtbdd_refs_pushptr(&cap_vars);
    for (int i=0; i<game.cap_count; i++) cap_vars = mtbdd_set_add(cap_vars, game.priobits+game.statebits+game.uap_count+i);

    // Construct ns_vars the set of next prio-state variables 
    MTBDD ns_vars = mtbdd_set_empty();
    mtbdd_refs_pushptr(&ns_vars);
    const int offset = game.cap_count + game.uap_count + game.priobits + game.statebits;
    for (int i=0; i<(game.priobits + game.statebits); i++) ns_vars = mtbdd_set_add(ns_vars, offset + i);

    // Compute the number of states from the strategy of the solved symbolic game (?)
    MTBDD states = sylvan_project(game.trans, s_vars);
    long noStates = (long)sylvan_satcount(states, s_vars);
    // std::cout << "there are " << noStates << " states" << std::endl;

    // mapper will store for each [decoded] state in the symbolic game, the parity game node id
    std::map<int, int> mapper;

    // Start constructing the parity game for given number of states
    pg::Game *pargame = new pg::Game(noStates);

    // First initialize a parity game vertex for every state in the symbolic game
    {
        int idx=0;

        uint8_t state_arr[game.statebits];
        MTBDD lf = mtbdd_enum_all_first(states, s_vars, state_arr, NULL);
        while (lf != mtbdd_false) {
            // decode state
            int state = 0;
            for (int i=0; i<game.statebits; i++) {
                state <<= 1;
                if (state_arr[i]) state |= 1;
            }

            pargame->init_vertex(idx, 0, 1, std::to_string(state)); // controlled by the Odd
            mapper[state] = idx++;

            lf = mtbdd_enum_all_next(states, s_vars, state_arr, NULL);
        }
    }

    // Check that the transition relation of the symbolic game does not have priobits on source states
    assert(mtbdd_getvar(game.trans) >= (unsigned) game.priobits);

    // index of next parity game vertex (for the priority on the edge)
    int vidx = noStates;

    std::vector<int> succ_state;  // for current state, the successors

    std::map<MTBDD, int> uap_to_vertex; // map after-uap to vertex
    std::map<MTBDD, int> cap_to_vertex; // map after-cap (priority) to vertex

    uint8_t state_arr[game.statebits];
    MTBDD lf = mtbdd_enum_all_first(game.trans, s_vars, state_arr, NULL);
    while (lf != mtbdd_false) {
        // decode state
        int _s = 0;
        for (int i=0; i<game.statebits; i++) {
            _s <<= 1;
            if (state_arr[i]) _s |= 1;
        }

        // translate _s to state
        int state = mapper.at(_s);

        // find all "successors" of this state after the environment plays (uap)
        std::set<MTBDD> after_uap;
        {
            uint8_t uap_arr[game.uap_count];
            MTBDD lf2 = mtbdd_enum_first(lf, uap_vars, uap_arr, NULL);
            while (lf2 != mtbdd_false) {
                after_uap.insert(lf2);
                lf2 = mtbdd_enum_next(lf, uap_vars, uap_arr, NULL);
            }
        }

        // start adding edges!
        pargame->e_start(state);

        for (MTBDD uap_bdd : after_uap) {
            // check if we have seen this symbolic state before
            auto search = uap_to_vertex.find(uap_bdd);
            if (search != uap_to_vertex.end()) {
                // already exists
                // add edge from <state> to <search->second>
                pargame->e_add(state, search->second);
            } else {
                int uapv = vidx++;
                pargame->init_vertex(uapv, 0, 0, std::to_string(uap_bdd)); // controlled by Even
                uap_to_vertex[uap_bdd] = uapv;
                // add edge from <state> to <uapv>
                pargame->e_add(state, uapv);
            }
        }

        // we're done adding edges!
        pargame->e_finish();

        lf = mtbdd_enum_all_next(game.trans, s_vars, state_arr, NULL);
    }

    // we now have all post-UAP in the uap_to_vertex map, so lets process them...

    for (auto x = uap_to_vertex.begin(); x != uap_to_vertex.end(); x++) {
        const MTBDD uap_bdd = x->first;
        const int uap_v = x->second;

        std::set<MTBDD> after_cap;
        {
            uint8_t cap_arr[game.cap_count];
            MTBDD lf = mtbdd_enum_first(uap_bdd, cap_vars, cap_arr, NULL);
            while (lf != mtbdd_false) {
                after_cap.insert(lf);
                lf = mtbdd_enum_next(uap_bdd, cap_vars, cap_arr, NULL);
            }
        }

        pargame->e_start(uap_v);

        for (MTBDD cap_bdd : after_cap) {
            // only if not yet there
            auto search2 = cap_to_vertex.find(cap_bdd);
            if (search2 != cap_to_vertex.end()) {
                pargame->e_add(uap_v, search2->second);
            } else {
                int capv = vidx++;

                uint8_t ns_arr[game.priobits+game.statebits];
                MTBDD lf2 = mtbdd_enum_all_first(cap_bdd, ns_vars, ns_arr, NULL);
                assert(lf2 == mtbdd_true);

                // decode priostate
                int pr = 0;
                for (int i=0; i<game.priobits; i++) {
                    pr <<= 1;
                    if (ns_arr[i]) pr |= 1;
                }

                pargame->init_vertex(capv, pr, 0, std::to_string(cap_bdd)); // controlled by any
                cap_to_vertex[cap_bdd] = capv;
                pargame->e_add(uap_v, capv);
            }
        }

        pargame->e_finish();
    }

    // we now have all post-CAP in the cap_to_vertex map, so lets process them...

    for (auto x = cap_to_vertex.begin(); x != cap_to_vertex.end(); x++) {
        const MTBDD cap_bdd = x->first;
        const int cap_v = x->second;

        uint8_t ns_arr[game.priobits+game.statebits];
        MTBDD lf2 = mtbdd_enum_all_first(cap_bdd, ns_vars, ns_arr, NULL);
        assert(lf2 == mtbdd_true);

        // decode target
        int to = 0;
        for (int i=0; i<game.statebits; i++) {
            to <<= 1;
            if (ns_arr[game.priobits+i]) to |= 1;
        }
        to = mapper.at(to);

        pargame->e_start(cap_v);
        pargame->e_add(cap_v, to);
        pargame->e_finish();
    }

    // tell Oink we're done adding stuff, resize game to final size
    pargame->v_resize(vidx);

    return pargame;
}


void
strategy_to_pg(SymGame *game)
{
    MTBDD s_vars = mtbdd_set_empty();
    mtbdd_refs_pushptr(&s_vars);
    for (int i=0; i<game->statebits; i++) s_vars = mtbdd_set_add(s_vars, game->priobits+i);

    const int offset = game->cap_count + game->uap_count + game->priobits + game->statebits;

    MTBDD ns_vars = mtbdd_set_empty();
    mtbdd_refs_pushptr(&ns_vars);
    for (int i=0; i<(game->priobits + game->statebits); i++) ns_vars = mtbdd_set_add(ns_vars, offset + i);

    MTBDD states = sylvan_project(game->strategies, s_vars);
    long noStates = (long)sylvan_satcount(states, s_vars);
    // std::cout << "there are " << noStates << " states" << std::endl;

    pg::Game *pargame = new pg::Game(noStates);

    std::map<int, int> mapper;
    {
        int idx=0;

        uint8_t state_arr[game->statebits];
        MTBDD lf = mtbdd_enum_all_first(states, s_vars, state_arr, NULL);
        while (lf != mtbdd_false) {
            // decode state
            int state = 0;
            for (int i=0; i<game->statebits; i++) {
                state <<= 1;
                if (state_arr[i]) state |= 1;
            }

            pargame->init_vertex(idx, 0, 1, std::to_string(state)); // controlled by the Odd
            mapper[state] = idx++;

            lf = mtbdd_enum_all_next(states, s_vars, state_arr, NULL);
        }
    }

    // DO A THING
    MTBDD full = sylvan_and(game->trans, game->strategies);
    // We only care about priority 0 priostates
    while (!mtbdd_isleaf(full) && mtbdd_getvar(full) < (unsigned) game->priobits) {
        full = mtbdd_getlow(full);
    }

    int vidx = noStates;

    std::vector<int> succ_state;  // for current state, the successors

    uint8_t state_arr[game->statebits];
    MTBDD lf = mtbdd_enum_all_first(full, s_vars, state_arr, NULL);
    while (lf != mtbdd_false) {
        // decode state
        int _s = 0;
        for (int i=0; i<game->statebits; i++) {
            _s <<= 1;
            if (state_arr[i]) _s |= 1;
        }

        int state = mapper.at(_s);

        // std::cout << "strategies for state " << state << ": " << lf << std::endl;

        MTBDD ns = sylvan_project(lf, ns_vars);
        mtbdd_refs_pushptr(&ns);
        
        uint8_t ns_arr[game->priobits+game->statebits];
        MTBDD lf2 = mtbdd_enum_all_first(ns, ns_vars, ns_arr, NULL);
        while (lf2 != mtbdd_false) {
            assert(lf2 == mtbdd_true);

            // decode priostate
            int pr = 0;
            for (int i=0; i<game->priobits; i++) {
                pr <<= 1;
                if (ns_arr[i]) pr |= 1;
            }
            int to = 0;
            for (int i=0; i<game->statebits; i++) {
                to <<= 1;
                if (ns_arr[game->priobits+i]) to |= 1;
            }
            to = mapper.at(to);

            pargame->init_vertex(vidx, pr, 1);
            pargame->e_start(vidx);
            pargame->e_add(vidx, to);
            pargame->e_finish();
            succ_state.push_back(vidx);
            vidx++;

            // std::cout << "to (" << pr << ") " << to << std::endl;
            
            lf2 = mtbdd_enum_all_next(ns, ns_vars, ns_arr, NULL);
        }

        pargame->e_start(state);
        for (int v : succ_state) pargame->e_add(state, v);
        pargame->e_finish();
        succ_state.clear();

        lf = mtbdd_enum_all_next(full, s_vars, state_arr, NULL);
    }


    // tell Oink we're done adding stuff, resize game to final size
    pargame->v_resize(vidx);
    pargame->write_pgsolver(std::cout);
}


class AIGmaker {
private:
    aiger *a;
    HoaData *data;
    SymGame *game;
    
    int lit; // current next literal

    int* uap_to_lit; // the input literal for each uncontrolled AP
    int* state_to_lit; // the latch literal for each state bit
    char** caps; // labels for controlled APs
    int* var_to_lit; // translate BDD variable (uap/state) to AIGER literal

    std::map<MTBDD, int> mapping; // map MTBDD to AIGER literal
    std::map<uint64_t, int> cache; // cache for ands

    int bdd_to_aig(MTBDD bdd);
    int makeand(int rhs0, int rhs1);

public:
    AIGmaker(HoaData *data, SymGame *game);
    ~AIGmaker();

    void processCAP(int i, MTBDD bdd);
    void processState(int i, MTBDD bdd);
    void write(FILE* out);
};


AIGmaker::AIGmaker(HoaData *data, SymGame *game) : data(data), game(game) {
    a = aiger_init();
    lit = 2;
    uap_to_lit = new int[game->uap_count];
    state_to_lit = new int[game->statebits];
    caps = new char*[game->cap_count];
    var_to_lit = new int[game->statebits+game->priobits+game->uap_count];

    // Set which APs are controllable in the bitset controllable
    pg::bitset controllable(data->noAPs);
    for (int i=0; i<data->noCntAPs; i++) {
        controllable[data->cntAPs[i]] = 1;
    }

    int uap_idx = 0;
    int cap_idx = 0;
    for (int i=0; i<data->noAPs; i++) {
        if (!controllable[i]) {
            uap_to_lit[uap_idx] = lit;
            var_to_lit[game->priobits+game->statebits+uap_idx] = lit;
            aiger_add_input(a, lit, data->aps[i]);
            uap_idx++;
            lit += 2;
        } else {
            caps[cap_idx++] = data->aps[i];
        }
    }

    for (int i=0; i<game->priobits; i++) {
        var_to_lit[i] = 0;
    }

    for (int i=0; i<game->statebits; i++) {
        state_to_lit[i] = lit;
        var_to_lit[game->priobits+i] = lit;
        lit += 2;
    }
}

AIGmaker::~AIGmaker()
{
    delete[] uap_to_lit;
    delete[] state_to_lit;
    delete[] caps;
    delete[] var_to_lit;
}


int
AIGmaker::makeand(int rhs0, int rhs1)
{
    if (rhs1 < rhs0) {
        int tmp = rhs0;
        rhs0 = rhs1;
        rhs1 = tmp;
    }

    if (rhs0 == 0) return 0;
    if (rhs0 == 1) return rhs1;

    uint64_t cache_key = rhs1;
    cache_key <<= 32;
    cache_key |= rhs0;
    auto c = cache.find(cache_key);
    if (c != cache.end()) {
        return c->second;
    } else {
        aiger_add_and(a, lit, rhs0, rhs1);
        cache[cache_key] = lit;
        lit += 2;
        return lit-2;
    }
}



int
AIGmaker::bdd_to_aig(MTBDD bdd)
{
    if (bdd == mtbdd_true) return aiger_true;
    if (bdd == mtbdd_false) return aiger_false;
 
    bool comp = false;
    if (bdd & sylvan_complement) {
        bdd ^= sylvan_complement;
        comp = true;
    }

    auto it = mapping.find(bdd);
    if (it != mapping.end()) {
        return comp ? aiger_not(it->second) : it->second;
    }

    int the_lit = var_to_lit[mtbdd_getvar(bdd)];

    MTBDD low = mtbdd_getlow(bdd);
    MTBDD high = mtbdd_gethigh(bdd);

    int res;

    if (low == mtbdd_false) {
        // only high (value 1)
        if (high == mtbdd_true) {
            // actually this is the end, just the lit
            res = the_lit;
        } else {
            // AND(the_lit, ...)
            int rhs0 = the_lit;
            int rhs1 = bdd_to_aig(high);
            res = makeand(rhs0, rhs1);
        }
    } else if (high == mtbdd_false) {
        // only low (value 0)
        if (low == mtbdd_true) {
            // actually this is the end, just the lit, negated
            res = aiger_not(the_lit);
        } else {
            // AND(not the_lit, ...)
            int rhs0 = aiger_not(the_lit);
            int rhs1 = bdd_to_aig(low);
            res = makeand(rhs0, rhs1);
        }
    } else {
        // OR(low, high) == ~AND(~AND(the_lit, ...), ~AND(~the_lit, ...))
        int lowres = bdd_to_aig(low);
        int highres = bdd_to_aig(high);
        int rhs0 = aiger_not(makeand(aiger_not(the_lit), lowres));
        int rhs1 = aiger_not(makeand(the_lit, highres));
        res = aiger_not(makeand(rhs0, rhs1));
    }
        
    mapping[bdd] = res;

    return comp ? aiger_not(res) : res;
}


void
AIGmaker::processCAP(int i, MTBDD bdd)
{
    int res = bdd_to_aig(bdd);
    aiger_add_output(a, res, caps[i]); // simple, really
}

void
AIGmaker::processState(int i, MTBDD bdd)
{
    int res = bdd_to_aig(bdd);
    aiger_add_latch(a, state_to_lit[i], res, "");
}

void
AIGmaker::write(FILE* out)
{
    aiger_write_to_file(a, aiger_ascii_mode, out);
}



bool
solveSymGame(SymGame *game)
{
    const int offset = game->cap_count + game->uap_count + game->priobits + game->statebits;

    MTBDD s_vars = mtbdd_set_empty();
    mtbdd_refs_pushptr(&s_vars);
    for (int i=0; i<(game->priobits + game->statebits); i++) s_vars = mtbdd_set_add(s_vars, i);

    MTBDD uap_vars = mtbdd_set_empty();
    mtbdd_refs_pushptr(&uap_vars);
    for (int i=0; i<game->uap_count; i++) uap_vars = mtbdd_set_add(uap_vars, game->priobits+game->statebits+i);

    MTBDD cap_vars = mtbdd_set_empty();
    mtbdd_refs_pushptr(&cap_vars);
    for (int i=0; i<game->cap_count; i++) cap_vars = mtbdd_set_add(cap_vars, game->priobits+game->statebits+game->uap_count+i);

    MTBDD ns_vars = mtbdd_set_empty();
    mtbdd_refs_pushptr(&ns_vars);
    for (int i=0; i<(game->priobits + game->statebits); i++) ns_vars = mtbdd_set_add(ns_vars, offset + i);

    /*
    MTBDD str_vars = mtbdd_set_empty();
    mtbdd_refs_pushptr(&str_vars);
    for (int i=0; i<game->statebits; i++) str_vars = mtbdd_set_add(str_vars, game->priobits + i);
    for (int i=0; i<game->uap_count; i++) str_vars = mtbdd_set_add(str_vars, game->priobits+game->statebits+i);
    for (int i=0; i<game->cap_count; i++) str_vars = mtbdd_set_add(str_vars, game->priobits+game->statebits+game->uap_count+i);
    for (int i=0; i<game->statebits; i++) str_vars = mtbdd_set_add(str_vars, offset + game->priobits + i);
    */

    /*
    MTBDD all_vars = mtbdd_set_empty();
    mtbdd_refs_pushptr(&all_vars);
    for (int i=0; i<(game->uap_count+game->cap_count+2*(game->priobits+game->statebits)); i++) all_vars = mtbdd_set_add(all_vars, i);
    */

    MTBDD odd = sylvan_ithvar(game->priobits-1); // deepest bit is the parity
    mtbdd_refs_pushptr(&odd);

    MTBDD from_next = mtbdd_map_empty();
    mtbdd_refs_pushptr(&from_next);
    for (int i=0; i<(game->priobits + game->statebits); i++) {
        from_next = mtbdd_map_add(from_next, offset+i, sylvan_ithvar(i));
    }

    MTBDD to_next = mtbdd_map_empty();
    mtbdd_refs_pushptr(&to_next);
    for (int i=0; i<(game->priobits + game->statebits); i++) {
        to_next = mtbdd_map_add(to_next, i, sylvan_ithvar(offset+i));
    }

    MTBDD* priostates = new MTBDD[game->maxprio+1];
    for (int i=0; i<=game->maxprio; i++) {
        priostates[i] = encodeprio(i, game->priobits);
        mtbdd_refs_pushptr(&priostates[i]);
    }

    MTBDD* lowereq = new MTBDD[game->maxprio+1];
    for (int i=0; i<=game->maxprio; i++) {
        lowereq[i] = mtbdd_false;
        mtbdd_refs_pushptr(&lowereq[i]);
        for (int j=0; j<=i; j++) {
            lowereq[i] = sylvan_or(priostates[j], lowereq[i]);
        }
    }

    // We force initial state to be state 0 in construction
    MTBDD initial = encode_priostate(0, 0, game->statebits, game->priobits, 0);
    mtbdd_refs_pushptr(&initial);

    // using the freezing FPI algorithm
    // 
    MTBDD distractions = mtbdd_false; // initially, no distractions
    mtbdd_refs_pushptr(&distractions);

    MTBDD strategies = mtbdd_false; // will have strategies ([prio]state -> priostate)
    mtbdd_refs_pushptr(&strategies);

    MTBDD* freeze = new MTBDD[game->maxprio+1];
    for (int i=0; i<=game->maxprio; i++) {
        freeze[i] = mtbdd_false;
        mtbdd_refs_pushptr(&freeze[i]);
    }

    int pr = 0;
    while (pr <= game->maxprio) {
        MTBDD onestepeven = mtbdd_false;
        mtbdd_refs_pushptr(&onestepeven); // push onestepeven

        {
            // Compute OneStepEven := (Distraction /\ OddPrio) \/ ~(Distraction \/ OddPrio)
            MTBDD OddAndDistraction = sylvan_and(distractions, odd);
            mtbdd_refs_pushptr(&OddAndDistraction);

            MTBDD OddOrDistraction = sylvan_or(distractions, odd);
            mtbdd_refs_pushptr(&OddOrDistraction);

            onestepeven = sylvan_or(OddAndDistraction, sylvan_not(OddOrDistraction));
            mtbdd_refs_popptr(2); // pop OddAndDistraction, OddOrDistraction

            // and convert to prime variables
            onestepeven = sylvan_compose(onestepeven, to_next);
        }

        // then take product with transition
        // remember the strat: state -> uap -> cap
        MTBDD strat = onestepeven = sylvan_and_exists(game->trans, onestepeven, ns_vars);
        mtbdd_refs_pushptr(&strat);

        // Extract vertices won by even in one step = \forall U. \exists C. onestepeven
        onestepeven = sylvan_exists(onestepeven, cap_vars);
        onestepeven = sylvan_forall(onestepeven, uap_vars);

        // Restrict strat to states that are one step even (winning)
        strat = sylvan_and(strat, onestepeven);

        // those are won by Even in one step ... 
        MTBDD newd = mtbdd_false;
        mtbdd_refs_pushptr(&newd); // push newd

        if ((pr&1) == 0) {
            // new EVEN distractions are states that are NOT one step even
            newd = sylvan_and(priostates[pr], sylvan_not(onestepeven));
        } else {
            // new ODD distractions are states that ARE one step even
            newd = sylvan_and(priostates[pr], onestepeven);
        }
        // remove vertices from newd that are already distractions
        newd = sylvan_and(newd, sylvan_not(distractions));

        if (newd != sylvan_false) {
            // add to distractions
            distractions = sylvan_or(distractions, newd);

            // now freeze/thaw/reset vertices <= pr
            MTBDD lwr = lowereq[pr];
            mtbdd_refs_pushptr(&lwr); // push lwr

            // remove vertices frozen at a higher priority
            for (int i=pr+1; i<=game->maxprio; i++) {
                // remove each higher freeze set from <lwr>
                lwr = sylvan_and(lwr, sylvan_not(freeze[i]));
            }

            // HANDLING THE STRATEGY
            //   strat := good_state -> uap -> cap
            // (only good/useful strategies) this way, we only store useful strategies
            // first restrict strat to <lwr> (<=pr and not frozen higher)
            strat = sylvan_and(strat, lwr);
            // next restrict strat to not frozen lower (so we preserve the correct strategy)
            for (int i=0; i<=pr; i++) {
                strat = sylvan_and(strat, sylvan_not(freeze[i]));
            }
            // now strat is only for unfrozen vertices in the <=pr game

            // Extract strategy (for controller) -- str_vars is all except priorities
            // strat = sylvan_project(strat, str_vars);
            MTBDD states_in_strat = sylvan_project(strat, s_vars);
            mtbdd_refs_pushptr(&states_in_strat);
            strategies = sylvan_ite(states_in_strat, strat, strategies); // big updater...
            mtbdd_refs_popptr(1); // pop states_in_strat

            // NOW: freeze!
            // same parity, also distraction = freeze
            if (pr&1) {
                // pr is odd, so freeze stuff won by even
                // keep lower odd distractions
                MTBDD oddlwr = sylvan_and(lwr, odd);
                mtbdd_refs_pushptr(&oddlwr);
                oddlwr = sylvan_and(oddlwr, distractions);
                // keep lower even nondistractions
                MTBDD evenlwr = sylvan_and(lwr, sylvan_not(odd));
                mtbdd_refs_pushptr(&evenlwr);
                evenlwr = sylvan_and(evenlwr, sylvan_not(distractions));
                // take union as new freeze set
                freeze[pr] = sylvan_or(oddlwr, evenlwr);
                mtbdd_refs_popptr(2); // pop oddlwr, evenlwr
            } else {
                // pr is even, so freeze stuff won by odd
                // keep lower odd nondistractions
                MTBDD oddlwr = sylvan_and(lwr, odd);
                mtbdd_refs_pushptr(&oddlwr);
                oddlwr = sylvan_and(oddlwr, sylvan_not(distractions));
                // keep lower even distractions
                MTBDD evenlwr = sylvan_and(lwr, sylvan_not(odd));
                mtbdd_refs_pushptr(&evenlwr);
                evenlwr = sylvan_and(evenlwr, distractions);
                // take union as new freeze set
                freeze[pr] = sylvan_or(oddlwr, evenlwr);
                mtbdd_refs_popptr(2); // pop oddlwr, evenlwr
            }

            // select all lower not frozen at pr and reset them
            lwr = sylvan_and(lwr, sylvan_not(freeze[pr])); 
            distractions = sylvan_and(distractions, sylvan_not(lwr));
            mtbdd_refs_popptr(1); // pop lwr

            // erase all lower freeze sets
            for (int i=0; i<pr; i++) freeze[i] = mtbdd_false;

            // finally reset pr to 0
            pr = 0;
        } else {
            // No distractions, so just update the strategy

            // Reduce strategy to <=pr game
            strat = sylvan_and(strat, lowereq[pr]);
            // Remove all frozen vertices from strat
            for (int i=0; i<=game->maxprio; i++) {
                strat = sylvan_and(strat, sylvan_not(freeze[i]));
            }
            // Now strat is only for unfrozen vertices in the <=pr game
            // Update <strategies> with <strat> but only for states in <strat>
            MTBDD states_in_strat = sylvan_project(strat, s_vars);
            mtbdd_refs_pushptr(&states_in_strat);
            strategies = sylvan_ite(states_in_strat, strat, strategies); // big updater...
            mtbdd_refs_popptr(1); // pop states_in_strat

            pr++;
        }

        mtbdd_refs_popptr(3); // pop newd, strat, onestepeven
    }

    // we now know if the initial state is distracting or not
    if (sylvan_and(initial, distractions) != sylvan_false) {
        mtbdd_refs_popptr(10+3*game->maxprio); // free it up

        delete[] freeze;
        delete[] priostates;
        delete[] lowereq;

        return false;
    }

    // We only care about priority 0 priostates
    while (!mtbdd_isleaf(strategies) && mtbdd_getvar(strategies) < (unsigned) game->priobits) {
        strategies = mtbdd_getlow(strategies);
    }

    // Select lowest strategy [heuristic]
    strategies = RUN(clarify, strategies, cap_vars);

    // Now remove all unreachable states according to the strategy  (slightly smaller controller)
    {
        MTBDD vars = mtbdd_set_empty();
        mtbdd_refs_pushptr(&vars);
        for (int i=0; i<game->priobits; i++) vars = mtbdd_set_add(vars, i);
        for (int i=0; i<game->priobits; i++) vars = mtbdd_set_add(vars, offset+i);
        vars = mtbdd_set_addall(vars, cap_vars);
        vars = mtbdd_set_addall(vars, uap_vars);

        MTBDD T = mtbdd_and_exists(strategies, game->trans, vars);
        mtbdd_refs_pushptr(&T);

        MTBDD visited = initial;
        mtbdd_refs_pushptr(&visited);

        MTBDD old = mtbdd_false;
        mtbdd_refs_pushptr(&old);

        // Just quick and dirty symbolic reachability ... this is not super efficient but it's OK
        while (old != visited) {
            old = visited;
            MTBDD next = mtbdd_and_exists(visited, T, s_vars); // cross product
            mtbdd_refs_pushptr(&next);
            next = sylvan_compose(next, from_next);   // and rename
            visited = sylvan_or(visited, next);
            mtbdd_refs_popptr(1);
        }

        strategies = sylvan_and(strategies, visited);
        mtbdd_refs_popptr(4);
    }

    game->strategies = strategies;

    {
        // compute bdds for the controllable APs
        for (int i=0; i<game->cap_count; i++) {
            MTBDD cap = sylvan_ithvar(game->uap_count+game->priobits+game->statebits+i);
            mtbdd_refs_pushptr(&cap);
            game->cap_bdds[i] = sylvan_and_exists(strategies, cap, cap_vars);
            mtbdd_refs_popptr(1);
        }
    }

    {
        // compute state bdds
        MTBDD su_vars = mtbdd_set_empty();
        mtbdd_refs_pushptr(&su_vars);
        for (int i=0; i<game->priobits; i++) su_vars = mtbdd_set_add(su_vars, i);
        for (int i=0; i<game->priobits; i++) su_vars = mtbdd_set_add(su_vars, offset+i);
        su_vars = mtbdd_set_addall(su_vars, cap_vars);

        MTBDD full = mtbdd_and_exists(strategies, game->trans, su_vars);
        mtbdd_refs_pushptr(&full);

        for (int i=0; i<game->statebits; i++) {
            MTBDD ns = sylvan_ithvar(offset+game->priobits+i);
            mtbdd_refs_pushptr(&ns);
            game->state_bdds[i] = sylvan_and_exists(full, ns, ns_vars);
            mtbdd_refs_popptr(1);
        }

        mtbdd_refs_popptr(2); // full and su_vars
    }


    mtbdd_refs_popptr(10+3*game->maxprio); // WHATEVER

    delete[] freeze;
    delete[] priostates;
    delete[] lowereq;

    return true;
}



cxxopts::ParseResult
handleOptions(int &argc, char**& argv)
{
    try {
        cxxopts::Options opts(argv[0], "HOA synthesis using Sylvan and Oink");
        opts.custom_help("[OPTIONS...] [FILE]");
        opts.add_options()
            ("help", "Print help")
            ("sym", "Generate and solve a symbolic parity game")
            ("naive", "Use the naive splitting procedure")
            ("print-game", "Just print the parity game")
            ("print-witness", "Print the witness parity game")
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
    const double t_before_parsing = wctime();
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

    const double t_after_parsing = wctime();
    if (verbose) {
        std::cerr << "finished parsing automaton in " << std::fixed << (t_after_parsing - t_before_parsing) << " sec." << std::endl;
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

    // Initialize Lace
    lace_start(1, 0); // initialize Lace, but sequentially

    // And initialize Sylvan
    sylvan_set_limits(128LL << 20, 1, 16); // should be enough (128 megabytes)
    sylvan_init_package();
    sylvan_init_mtbdd();

    bool explicit_solver = options["sym"].count() == 0;
    bool naive_splitting = options["naive"].count() > 0;
    bool write_pg = options["print-game"].count() > 0;

    if (explicit_solver) {
        const double t_before_splitting = wctime();
        // Remember the start vertex
        int vstart = data->start[0];

        // Construct the game
        pg::Game *game;
        if (naive_splitting) {
            game = constructGameNaive(data, isMaxParity, controllerIsOdd);
        } else {
            game = constructGame(data, isMaxParity, controllerIsOdd);
        }
        game->set_label(vstart, "initial");
        const double t_after_splitting = wctime();

        // free HOA allocated data structure
        resetHoa(data);
        if (verbose) std::cerr << "finished constructing game in " << std::fixed << (t_after_splitting - t_before_splitting) << " sec." << std::endl;

        if (verbose) sylvan_stats_report(stdout);

        // We don't need Sylvan anymore at this point
        sylvan_quit();

        if (write_pg) {
            // in case we want to write the file to PGsolver file format...
            game->write_pgsolver(std::cout);
            // std::cerr << "initial vertex: " << vstart << std::endl;
            exit(0);
        }

        if (verbose) {
            std::cerr << "constructed game has " << game->vertexcount() << " vertices and " << game->edgecount() << " edges." << std::endl;
        }

        // we sort now, so we can track the initial state
        double begin = wctime();
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
        engine.run();
        double end = wctime();

        // report how long it all took
        if (verbose) std::cerr << "finished solving game in " << std::fixed << (end-begin) << " sec." << std::endl;

        // finally, check if the initial vertex is won by controller or environment
        if (game->winner[vstart] == 0) {
            std::cout << "REALIZABLE";
            exit(10);
        } else {
            std::cout << "UNREALIZABLE";
            exit(20);
        }
    } else {
        // Construct the game
        const double t_before_construct = wctime();
        auto sym = constructSymGame(data, isMaxParity, controllerIsOdd);
        const double t_after_construct = wctime();

        if (write_pg) {
            // in case we want to write the file to PGsolver file format...
            auto pg = symgame_to_pg(*sym);
            pg->write_pgsolver(std::cout);
            exit(0);
        }

        if (verbose) std::cerr << "finished constructing game in " << std::fixed << (t_after_construct - t_before_construct) << " sec." << std::endl;

        const double t_before_solve = wctime();
        bool res = solveSymGame(sym);
        const double t_after_solve = wctime();

        if (verbose) std::cerr << "finished solving game in " << std::fixed << (t_after_solve - t_before_solve) << " sec." << std::endl;

        if (verbose) sylvan_stats_report(stdout);

        if (res) {
            /** 
             * maybe print witness parity game, which should be fully won by even
             */

            if (options["print-witness"].count() > 0) {
                strategy_to_pg(sym);
                exit(10);
            }

            AIGmaker maker(data, sym);
            for (int i=0; i<sym->cap_count; i++) {
                maker.processCAP(i, sym->cap_bdds[i]);
            }
            for (int i=0; i<sym->statebits; i++) {
                maker.processState(i, sym->state_bdds[i]);
            }

            std::cout << "REALIZABLE" << std::endl;
            maker.write(stdout);
            exit(10);
        } else {
            std::cout << "UNREALIZABLE" << std::endl;
            exit(20);
        }

        // free HOA allocated data structure
        resetHoa(data);
    }

    lace_stop();
}

