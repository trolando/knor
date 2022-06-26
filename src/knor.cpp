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
#include <deque>

#include <game.hpp>
#include <oink.hpp>
#include <solvers.hpp>
#include <tools/cxxopts.hpp>
#include <sylvan.h>
#include <knor.hpp>
#include <symgame.hpp>
#include <bisim.hpp>

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
encode_state(uint32_t state, const int statebits, const int priobits, const int offset)
{
    // create a cube
    MTBDD cube = mtbdd_true;
    for (int i=0; i<statebits; i++) {
        const int bit = offset+statebits+priobits-i-1;
        if (state & 1) cube = mtbdd_makenode(bit, mtbdd_false, cube);
        else cube = mtbdd_makenode(bit, cube, mtbdd_false);
        state >>= 1;
    }
    return cube;
}


/**
* Encode a priostate as a BDD, using priobits > statebits
* High-significant bits come before low-significant bits in the BDD
*/
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
collect_targets(MTBDD trans, std::set<uint64_t> &res, const int statebits, const int priobits)
{
    if (mtbdd_isleaf(trans)) {
        uint64_t leaf = mtbdd_getint64(trans);
        res.insert(leaf);

        uint32_t priority = (uint32_t)(leaf>>32);
        uint32_t state = (uint32_t)(leaf & 0xffffffff);

        return encode_priostate(state, priority, statebits, priobits, 0);
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
            lblbdd = RUN(evalLabel, label, data, variables);
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


class AIGmaker {
private:
    aiger *a;
    HoaData *data;
    SymGame *game;

    bool isop = false; // use ISOP
    bool verbose = false;
    
    int lit; // current next literal

    int* uap_to_lit; // the input literal for each uncontrolled AP
    int* state_to_lit; // the latch literal for each state bit
    char** caps; // labels for controlled APs
    std::map<uint32_t, int> var_to_lit; // translate BDD variable (uap/state) to AIGER literal

    MTBDD* cap_bdds;   // contains the solution: controllable ap bdds: state -> uap -> B
    MTBDD* state_bdds; // contains the solution: state bit bdds      : state -> uap -> B

    std::map<MTBDD, int> mapping; // map MTBDD to AIGER literal
    std::map<uint64_t, int> cache; // cache for ands

    int makeand(int rhs0, int rhs1);
    int bdd_to_aig(MTBDD bdd);           // use recursive encoding of BDD (shannon expanion)
    int bdd_to_aig_isop(MTBDD bdd);
    int bdd_to_aig_cover(ZDD bdd);       // use recursive encoding of ZDD cover (~shannon expansion)
    int bdd_to_aig_cover_sop(ZDD cover); // use SOP encoding ("two level logic")
    void simplify_and(std::deque<int> &gates);
    void simplify_or(std::deque<int> &gates);

    void processCAP(int i, MTBDD bdd);
    void processState(int i, MTBDD bdd);

public:
    AIGmaker(HoaData *data, SymGame *game);
    ~AIGmaker();

    void setIsop()
    {
        this->isop = true;
    }

    void setVerbose()
    {
        this->verbose = true;
    }

    void process();
    void write(FILE* out);
};


AIGmaker::AIGmaker(HoaData *data, SymGame *game) : data(data), game(game)
{
    a = aiger_init();
    lit = 2;
    uap_to_lit = new int[game->uap_count];
    state_to_lit = new int[game->statebits];
    caps = new char*[game->cap_count];

    // Set which APs are controllable in the bitset controllable
    pg::bitset controllable(data->noAPs);
    for (int i=0; i<data->noCntAPs; i++) {
        controllable[data->cntAPs[i]] = 1;
    }

    cap_bdds = new MTBDD[game->cap_count];
    for (int i=0; i<game->cap_count; i++) {
        cap_bdds[i] = mtbdd_false;
        mtbdd_protect(&cap_bdds[i]);
    }

    state_bdds = new MTBDD[game->statebits];
    for (int i=0; i<game->statebits; i++) {
        state_bdds[i] = mtbdd_false;
        mtbdd_protect(&state_bdds[i]);
    }

    int uap_idx = 0;
    int cap_idx = 0;
    for (int i=0; i<data->noAPs; i++) {
        if (!controllable[i]) {
            uap_to_lit[uap_idx] = lit;
            aiger_add_input(a, lit, data->aps[i]);
            uap_idx++;
            lit += 2;
        } else {
            caps[cap_idx++] = data->aps[i];
        }
    }

    for (int i=0; i<game->statebits; i++) {
        state_to_lit[i] = lit;
        lit += 2;
    }

    {
        MTBDD _vars = game->s_vars;
        for (int i=0; i<game->statebits; i++) {
            uint32_t bddvar = mtbdd_set_first(_vars);
            _vars = mtbdd_set_next(_vars);
            var_to_lit[bddvar] = state_to_lit[i];
            // std::cerr << "BDD variable " << bddvar << " is state " << i << " literal " << state_to_lit[i] << std::endl;
        }
        _vars = game->uap_vars;
        for (int i=0; i<game->uap_count; i++) {
            uint32_t bddvar = mtbdd_set_first(_vars);
            _vars = mtbdd_set_next(_vars);
            var_to_lit[bddvar] = uap_to_lit[i];
            // std::cerr << "BDD variable " << bddvar << " is uAP " << i << " literal " << uap_to_lit[i] << std::endl;
        }
    }

    {
        // compute bdds for the controllable APs
        for (int i=0; i<game->cap_count; i++) {
            MTBDD cap = sylvan_ithvar(mtbdd_set_first(game->cap_vars)+i);
            mtbdd_protect(&cap);
            // keep just s and u... get rid of other cap variables
            this->cap_bdds[i] = sylvan_and_exists(game->strategies, cap, game->cap_vars);
            mtbdd_unprotect(&cap);
        }
    }

    {
        // compute state bdds
        MTBDD su_vars = mtbdd_set_addall(game->p_vars, game->cap_vars);
        mtbdd_protect(&su_vars);

        // su_vars is priority and cap, which is to be removed...
        // full is: s > u > ns
        MTBDD full = mtbdd_and_exists(game->strategies, game->trans, su_vars);
        mtbdd_protect(&full);

        MTBDD ns = mtbdd_false;
        mtbdd_protect(&ns);

        for (int i=0; i<game->statebits; i++) {
            ns = sylvan_ithvar(mtbdd_set_first(game->ns_vars) + i);
            // don't care about the other next state variables... keep just s and u
            this->state_bdds[i] = sylvan_and_exists(full, ns, game->ns_vars);
        }

        mtbdd_unprotect(&ns);
        mtbdd_unprotect(&full);
        mtbdd_unprotect(&su_vars);
    }
}

AIGmaker::~AIGmaker()
{
    delete[] uap_to_lit;
    delete[] state_to_lit;
    delete[] caps;

    for (int i=0; i<game->cap_count; i++) mtbdd_unprotect(&cap_bdds[i]);
    for (int i=0; i<game->statebits; i++) mtbdd_unprotect(&state_bdds[i]);

    delete[] cap_bdds;
    delete[] state_bdds;
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
    if (1 and c != cache.end()) {
        return c->second;
    } else {
        aiger_add_and(a, lit, rhs0, rhs1);
        cache[cache_key] = lit;
        lit += 2;
        return lit-2;
    }
}

void
AIGmaker::simplify_and(std::deque<int> &gates)
{
    // for each pair of gates in gates, check the cache
    for (auto first = gates.begin(); first != gates.end(); ++first) {
        for (auto second = first + 1; second != gates.end(); ++second) {
            int left = *first;
            int right = *second;
            if (left > right) std::swap(left, right);
            uint64_t cache_key = right;
            cache_key <<= 32;
            cache_key |= left;
            auto c = cache.find(cache_key);
            if (c != cache.end()) {
                //gates.erase(std::remove_if(gates.begin(), gates.end(), [=](int x){return x==left or x==right;}),
                //        gates.end());
                gates.erase(second);
                gates.erase(first);
                gates.push_back(c->second);
                simplify_and(gates);
                return;
            }
        }
    }
}

void
AIGmaker::simplify_or(std::deque<int> &gates)
{
    // for each pair of gates in gates, check the cache
    for (auto first = gates.begin(); first != gates.end(); ++first) {
        for (auto second = first + 1; second != gates.end(); ++second) {
            int left = aiger_not(*first);
            int right = aiger_not(*second);
            if (left > right) std::swap(left, right);
            uint64_t cache_key = right;
            cache_key <<= 32;
            cache_key |= left;
            auto c = cache.find(cache_key);
            if (c != cache.end()) {
                gates.erase(second);
                gates.erase(first);
                gates.push_back(aiger_not(c->second));
                simplify_or(gates);
                return;
            }
        }
    }
}

int
AIGmaker::bdd_to_aig_isop(MTBDD bdd)
{
    if (verbose) {
        std::cerr << "running isop for BDD with " << mtbdd_nodecount(bdd) << " nodes." << std::endl;
    }
    MTBDD bddres;
    ZDD isop = zdd_isop(bdd, bdd, &bddres);
    zdd_protect(&isop);
    // no need to reference the result...
    assert(bdd == bddres);
    assert(bdd == zdd_cover_to_bdd(isop));
    if (verbose) {
        std::cerr << "isop has " << (long)zdd_pathcount(isop) << " terms and " << zdd_nodecount(&isop, 1) << " nodes." << std::endl;
    }

    int res = bdd_to_aig_cover(isop);
    zdd_unprotect(&isop);
    return res;
}

int
AIGmaker::bdd_to_aig_cover_sop(ZDD cover)
{
    if (cover == zdd_true) return aiger_true;
    if (cover == zdd_false) return aiger_false;

    // a product could consist of all variables, and a -1 to denote
    //  the end of the product
    int product[game->statebits+game->uap_count+1] = { 0 };

    // a queue that stores all products, which will need to be summed
    std::deque<int> products;

    ZDD res = zdd_cover_enum_first(cover, product);
    while (res != zdd_false) {
        //  containing subproducts in the form of gates
        std::deque<int> gates;

        for (int i=0; product[i] != -1; i++) {
            int the_lit = var_to_lit[product[i]/2];
            if (product[i]&1) the_lit = aiger_not(the_lit);
            gates.push_back(the_lit);
        }

        // simplify_and(gates);

        // while we still have subproducts we need to AND together
        while (!gates.empty()) {
            int last = gates.front();
            gates.pop_front();
            if (!gates.empty()) {
                int last2 = gates.front();
                gates.pop_front();
                int new_gate = makeand(last, last2);
                gates.push_back(new_gate);
            } else {
                products.push_back(last);
            }
        }
        res = zdd_cover_enum_next(cover, product); // go to the next product
    }

    // products queue should now be full of complete products that need to be summed

    // simplify_or(products);

    while (!products.empty()) {
        int product1 = products.front();
        products.pop_front();
        if (!products.empty()) {
            int product2 = products.front();
            products.pop_front();
            int summed_product = aiger_not(makeand(aiger_not(product1), aiger_not(product2)));
            products.push_back(summed_product);
        } else { // product1 is the final sum of all products
            return product1;
        }
    }

    return aiger_false; // should be unreachable, tbh.
}

int
AIGmaker::bdd_to_aig_cover(ZDD cover)
{
    if (cover == zdd_true) return aiger_true;
    if (cover == zdd_false) return aiger_false;

    auto it = mapping.find(cover);
    if (it != mapping.end()) {
        return it->second;
    }

    int the_var = zdd_getvar(cover);
    int the_lit = var_to_lit[the_var/2];
    if (the_var & 1) the_lit = aiger_not(the_lit);

    ZDD low = zdd_getlow(cover);
    ZDD high = zdd_gethigh(cover);

    int res = the_lit;

    if (high != zdd_true) {
        auto x = bdd_to_aig_cover(high);
        res = makeand(res, x);
    }

    if (low != zdd_false) {
        auto x = bdd_to_aig_cover(low);
        res = aiger_not(makeand(aiger_not(res), aiger_not(x)));
    }

    mapping[cover] = res;
    return res;
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
    int res = isop ? bdd_to_aig_isop(bdd) : bdd_to_aig(bdd);
    aiger_add_output(a, res, caps[i]); // simple, really
}

void
AIGmaker::processState(int i, MTBDD bdd)
{
    int res = isop ? bdd_to_aig_isop(bdd) : bdd_to_aig(bdd);
    aiger_add_latch(a, state_to_lit[i], res, "");
}

void
AIGmaker::process()
{
    // if ISOP, first convert all cap bdds etc to covers
    if (isop) {
        ZDD* cap_zdds = new ZDD[game->cap_count];
        for (int i=0; i<game->cap_count; i++) {
            cap_zdds[i] = zdd_false;
            zdd_protect(&cap_zdds[i]);

            MTBDD bddres;
            cap_zdds[i] = zdd_isop(cap_bdds[i], cap_bdds[i], &bddres);
            assert(bddres == cap_bdds[i]);
            if (verbose) {
                std::cerr << "isop has " << (long)zdd_pathcount(cap_zdds[i]) << " terms and " << zdd_nodecount(&cap_zdds[i], 1) << " nodes." << std::endl;
            }
        }

        ZDD* state_zdds = new ZDD[game->statebits];
        for (int i=0; i<game->statebits; i++) {
            state_zdds[i] = zdd_false;
            zdd_protect(&state_zdds[i]);

            MTBDD bddres;
            state_zdds[i] = zdd_isop(state_bdds[i], state_bdds[i], &bddres);
            assert(bddres == state_bdds[i]);
            if (verbose) {
                std::cerr << "isop has " << (long)zdd_pathcount(state_zdds[i]) << " terms and " << zdd_nodecount(&state_zdds[i], 1) << " nodes." << std::endl;
            }
        }

        for (int i=0; i<game->cap_count; i++) {
            int res = bdd_to_aig_cover(cap_zdds[i]);
            aiger_add_output(a, res, caps[i]); // simple, really
        }
        for (int i=0; i<game->statebits; i++) {
            int res = bdd_to_aig_cover(state_zdds[i]);
            aiger_add_latch(a, state_to_lit[i], res, "");
        }
    } else {
        for (int i=0; i<game->cap_count; i++) {
            processCAP(i, cap_bdds[i]);
        }
        for (int i=0; i<game->statebits; i++) {
            processState(i, state_bdds[i]);
        }
    }
}

void
AIGmaker::write(FILE* out)
{
    aiger_write_to_file(a, aiger_ascii_mode, out);
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
            ("explicit", "Use the explicit splitting procedure")
            ("isop", "Generate AIG using ISOP instead of ITE")
            ("bisim", "Apply bisimulation minimisation to the solution")
            ("real", "Only check realiziability, don't synthesize")
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


VOID_TASK_0(gc_start)
{
    std::cerr << "starting garbage collection..." << std::endl;
}


VOID_TASK_0(gc_end)
{
    std::cerr << "garbage collection finished." << std::endl;
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
    sylvan_set_limits(512LL << 20, 1, 14); // should be enough (512 megabytes)
    sylvan_init_package();
    sylvan_init_mtbdd();
    sylvan_init_zdd();
    if (verbose) sylvan_gc_hook_pregc(TASK(gc_start));
    if (verbose) sylvan_gc_hook_postgc(TASK(gc_end));

    bool explicit_solver = options["sym"].count() == 0;
    bool naive_splitting = options["naive"].count() > 0;
    bool explicit_splitting = options["explicit"].count() > 0;
    bool write_pg = options["print-game"].count() > 0;

    SymGame *sym = nullptr;
    bool realizable = false; // not known yet

    if (explicit_solver) {
        // Remember the start vertex
        int vstart = data->start[0];
        std::map<int, MTBDD> vertex_to_bdd;

        // Construct the explicit game
        pg::Game *game = nullptr;
        const double t_before_splitting = wctime();
        if (naive_splitting) {
            game = RUN(constructGameNaive, data, isMaxParity, controllerIsOdd);
        } else if (explicit_splitting) {
            game = RUN(constructGame, data, isMaxParity, controllerIsOdd);
        } else {
            vstart = 0; // always set to 0 by constructSymGame
            sym = RUN(constructSymGame, data, isMaxParity, controllerIsOdd);
            game = sym->toExplicit(vertex_to_bdd);
        }
        const double t_after_splitting = wctime();

        if (verbose) {
            std::cerr << "finished constructing game in " << std::fixed << (t_after_splitting - t_before_splitting) << " sec." << std::endl;
            std::cerr << "constructed game has " << game->vertexcount() << " vertices and " << game->edgecount() << " edges." << std::endl;
        }

        if (write_pg) {
            // in case we want to write the file to PGsolver file format...
            game->set_label(vstart, "initial");
            game->write_pgsolver(std::cout);
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
        engine.setTrace(verbose ? 1 : 0);
        engine.setRenumber();
        engine.setSolver(solver);
        engine.setWorkers(-1);

        // and run the solver
        engine.run();
        double end = wctime();

        // report how long it all took
        if (verbose) std::cerr << "finished solving game in " << std::fixed << (end-begin) << " sec." << std::endl;

        // undo the sorting
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
        sym = RUN(constructSymGame, data, isMaxParity, controllerIsOdd);
        const double t_after_construct = wctime();

        if (write_pg) {
            // in case we want to write the file to PGsolver file format...
            std::map<int, MTBDD> vertex_to_bdd;
            auto pg = sym->toExplicit(vertex_to_bdd);
            pg->write_pgsolver(std::cout);
            exit(0);
        }

        if (verbose) std::cerr << "finished constructing game in " << std::fixed << (t_after_construct - t_before_construct) << " sec." << std::endl;

        const double t_before_solve = wctime();
        realizable = RUN(wrap_solve, sym, verbose);
        const double t_after_solve = wctime();

        if (verbose) std::cerr << "finished solving game in " << std::fixed << (t_after_solve - t_before_solve) << " sec." << std::endl;
    }

    if (options["real"].count() > 0) {
        if (realizable) {
            std::cout << "REALIZABLE" << std::endl;
            exit(10);
        } else {
            std::cout << "UNREALIZABLE" << std::endl;
            exit(20);
        }
    }

    if (realizable) {
        if (verbose) {
            // sym->print_trans();
            // sym->print_strategies();
        }

        RUN(wrap_pp, sym, verbose);

        if (verbose) {
            // sym->print_trans();
            // sym->print_strategies();
        }

        if (options["bisim"].count() > 0) {
            MTBDD partition = RUN(min_lts_strong, sym);
            mtbdd_protect(&partition);

            if (verbose) {
                // RUN(print_partition, sym, partition);
            }

            RUN(minimize, sym, partition, verbose);
            mtbdd_unprotect(&partition);
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
            exit(10);
        }

        // sylvan_gc();

        AIGmaker maker(data, sym);
        if (verbose) {
            maker.setVerbose();
        }
        if (options["isop"].count() > 0) {
            maker.setIsop();
        }
        maker.process();
        std::cout << "REALIZABLE" << std::endl;
        maker.write(stdout);
        exit(10);
    } else {
        std::cout << "UNREALIZABLE" << std::endl;
        exit(20);
    }

    if (verbose) sylvan_stats_report(stdout);

    // We don't need Sylvan anymore at this point
    sylvan_quit();

    // free HOA allocated data structure
    resetHoa(data);

    lace_stop();
}

