/**************************************************************************
 * Copyright (c) 2020-2021 Tom van Dijk
 *************************************************************************/

#include <iostream>
#include <map>
#include <set>
#include <deque>

#include <oink/game.hpp>
#include <oink/solvers.hpp>
#include <sylvan.h>
#include <knor.hpp>
#include <symgame.hpp>
#include <explicitgame.hpp>

extern "C" {
#include "simplehoa.h"
#include "aiger.h"
}

using namespace sylvan;

SymGame::SymGame(int statebits, int priobits, int uap_count, int cap_count, int maxprio)
{
    trans = mtbdd_false;
    strategies = mtbdd_false;

    mtbdd_protect(&trans);
    mtbdd_protect(&strategies);

    p_vars = mtbdd_set_empty();
    s_vars = mtbdd_set_empty();
    uap_vars = mtbdd_set_empty();
    cap_vars = mtbdd_set_empty();
    np_vars = mtbdd_set_empty();
    ns_vars = mtbdd_set_empty();
    ps_vars = mtbdd_set_empty();
    pns_vars = mtbdd_set_empty();
    uns_vars = mtbdd_set_empty();

    mtbdd_protect(&p_vars);
    mtbdd_protect(&s_vars);
    mtbdd_protect(&uap_vars);
    mtbdd_protect(&cap_vars);
    mtbdd_protect(&np_vars);
    mtbdd_protect(&ns_vars);
    mtbdd_protect(&ps_vars);
    mtbdd_protect(&pns_vars);
    mtbdd_protect(&uns_vars);

    // skip to priobits
    int var = 0;
    for (int i=0; i<priobits; i++) p_vars = mtbdd_set_add(p_vars, var++);
    for (int i=0; i<statebits; i++) s_vars = mtbdd_set_add(s_vars, var++);
    for (int i=0; i<uap_count; i++) uap_vars = mtbdd_set_add(uap_vars, var++);
    for (int i=0; i<cap_count; i++) cap_vars = mtbdd_set_add(cap_vars, var++);
    for (int i=0; i<priobits; i++) np_vars = mtbdd_set_add(np_vars, var++);
    for (int i=0; i<statebits; i++) ns_vars = mtbdd_set_add(ns_vars, var++);
    ps_vars = mtbdd_set_addall(p_vars, s_vars);
    pns_vars = mtbdd_set_addall(np_vars, ns_vars);
    uns_vars = mtbdd_set_addall(uap_vars, ns_vars);

    this->maxprio = maxprio;
    this->priobits = priobits;
    this->statebits = statebits;
    this->uap_count = uap_count;
    this->cap_count = cap_count;
}

SymGame::~SymGame()
{
    mtbdd_unprotect(&trans);
    mtbdd_unprotect(&strategies);
    mtbdd_unprotect(&p_vars);
    mtbdd_unprotect(&s_vars);
    mtbdd_unprotect(&uap_vars);
    mtbdd_unprotect(&cap_vars);
    mtbdd_unprotect(&np_vars);
    mtbdd_unprotect(&ns_vars);
    mtbdd_unprotect(&ps_vars);
    mtbdd_unprotect(&pns_vars);
    mtbdd_unprotect(&uns_vars);
}

/**
 * Encode a state as a BDD, using statebits 0..<statebits>, offsetted by <offset>+<priobits>
 * High-significant bits come before low-significant bits in the BDD
 */
MTBDD SymGame::encode_state(uint32_t state, MTBDD statevars)
{
    // convert statevars to stack
    std::vector<unsigned int> vars;
    while (statevars != mtbdd_set_empty()) {
        vars.push_back(mtbdd_set_first(statevars));
        statevars = mtbdd_set_next(statevars);
    }
    // create a cube
    auto cube = mtbdd_true;
    while (!vars.empty()) {
        auto var = vars.back();
        vars.pop_back();
        if (state & 1) cube = mtbdd_makenode(var, mtbdd_false, cube);
        else cube = mtbdd_makenode(var, cube, mtbdd_false);
        state >>= 1;
    }
    return cube;
}

/**
 * Encode priority i.e. all states via priority <priority>
 */
MTBDD SymGame::encode_prio(int priority, int priobits)
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
* Encode a priostate as a BDD, with priobits before statebits
* High-significant bits come before low-significant bits in the BDD
*/
MTBDD SymGame::encode_priostate(uint32_t state, uint32_t priority, MTBDD statevars, MTBDD priovars)
{
    // convert statevars to stack
    std::vector<unsigned int> vars;
    while (statevars != mtbdd_set_empty()) {
        vars.push_back(mtbdd_set_first(statevars));
        statevars = mtbdd_set_next(statevars);
    }
    // create the cube
    auto cube = mtbdd_true;
    while (!vars.empty()) {
        auto var = vars.back();
        vars.pop_back();
        if (state & 1) cube = mtbdd_makenode(var, mtbdd_false, cube);
        else cube = mtbdd_makenode(var, cube, mtbdd_false);
        state >>= 1;
    }
    // convert priovars to stack
    while (priovars != mtbdd_set_empty()) {
        vars.push_back(mtbdd_set_first(priovars));
        priovars = mtbdd_set_next(priovars);
    }
    // create the rest of the cube
    while (!vars.empty()) {
        auto var = vars.back();
        vars.pop_back();
        if (priority & 1) cube = mtbdd_makenode(var, mtbdd_false, cube);
        else cube = mtbdd_makenode(var, cube, mtbdd_false);
        priority >>= 1;
    }
    return cube;
}

/**
 * Convert a transition label (Btree) to a BDD encoding the label
 * a label is essentially a boolean combination of atomic propositions and aliases
 */
TASK_3(MTBDD, evalLabel, BTree*, label, HoaData*, data, uint32_t*, variables)
{
    MTBDD left;
    MTBDD right;
    MTBDD result;
    switch (label->type) {
        case NT_BOOL:
            return label->id ? mtbdd_true : mtbdd_false;
        case NT_AND:
            left = CALL(evalLabel, label->left, data, variables);
            mtbdd_refs_pushptr(&left);
            right = CALL(evalLabel, label->right, data, variables);
            mtbdd_refs_pushptr(&right);
            result = sylvan_and(left, right);
            mtbdd_refs_popptr(2);
            return result;
        case NT_OR:
            left = CALL(evalLabel, label->left, data, variables);
            mtbdd_refs_pushptr(&left);
            right = CALL(evalLabel, label->right, data, variables);
            mtbdd_refs_pushptr(&right);
            result = sylvan_or(left, right);
            mtbdd_refs_popptr(2);
            return result;
        case NT_NOT:
            left = CALL(evalLabel, label->left, data, variables);
            mtbdd_refs_pushptr(&left);
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
 * Construct the symbolic game
 */
std::unique_ptr<SymGame> SymGame::constructSymGame(HoaData *data, bool isMaxParity, bool controllerIsOdd) {
    int vstart = data->start[0];

    // Set which APs are controllable in the bitset controllable
    pg::bitset controllable(data->noAPs);
    for (int i=0; i<data->noCntAPs; i++) {
        controllable[data->cntAPs[i]] = true;
    }

    // Count the number of controllable/uncontrollable APs
    const auto cap_count = controllable.count();
    const auto uap_count = data->noAPs - cap_count;

    // Prepare number of statebits and priobits
    int statebits = 1;
    while ((1ULL<<(statebits)) < (uint64_t)data->noStates) statebits++;

    const int evenMax = 2 + 2*((data->noAccSets+1)/2); // should be enough...
    int priobits = 1;
    while ((1ULL<<(priobits)) <= (unsigned)evenMax) priobits++;

    // Initially, set maxprio to 0, will be updated later
    auto res = std::make_unique<SymGame>(statebits, priobits, uap_count, cap_count, 0);

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
            BTree* label = state->label != nullptr ? state->label : trans->label;
            assert(label != nullptr);
            // there should be a priority at state or transition level
            int priority;
            if (trans->noAccSig != 0) {
                // there should be exactly one acceptance set!
                assert(trans->noAccSig == 1);
                // adjust priority
                priority = ExplicitGame::adjustPriority(trans->accSig[0], isMaxParity, controllerIsOdd, data->noAccSets);
            } else {
                auto id = trans->successors[0];
                assert(data->states[id].noAccSig == 1);
                priority = ExplicitGame::adjustPriority(data->states[id].accSig[0], isMaxParity, controllerIsOdd, data->noAccSets);
            }
            if (priority > res->maxprio) res->maxprio = priority;
            // swap with initial if needed
            int succ = trans->successors[0];
            if (succ == 0) succ = vstart;
            else if (succ == vstart) succ = 0;
            // encode the label as a MTBDD
            lblbdd = RUN(evalLabel, label, data, variables);
            // encode priostate (leaf) and update transition relation
            leaf = encode_priostate(succ, priority, res->ns_vars, res->np_vars);
            // trans := lbl THEN leaf ELSE trans
            transbdd = mtbdd_ite(lblbdd, leaf, transbdd);
            // deref lblbdd and leaf
            lblbdd = leaf = mtbdd_false;
        }

        // encode source state and add to full transition relation
        int src = state->id;
        if (src == 0) src = vstart;
        else if (src == vstart) src = 0;
        statebdd = encode_state(src, res->s_vars);
        // update full trans with <statebdd> then <transbdd>
        res->trans = mtbdd_ite(statebdd, transbdd, res->trans);
        // deref statebdd and transbdd
        statebdd = transbdd = mtbdd_false;
    }

    mtbdd_refs_popptr(4);

    return res;
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
            mtbdd_refs_push(low);
            MTBDD high = CALL(clarify, mtbdd_gethigh(str), cap_vars);
            mtbdd_refs_pop(1);
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


/**
 * Convert the symbolic parity game to an explicit parity game, labeling vertices with BDD node indices.
 */
std::unique_ptr<pg::Game> SymGame::toExplicit(std::map<int, MTBDD> &vertex_to_bdd)
{
    // Check that the transition relation of the symbolic game does not have priobits on source states
    assert(mtbdd_getvar(this->trans) >= (unsigned) mtbdd_set_first(s_vars));

    // Compute the number of states with transitions
    auto states = sylvan_project(this->trans, s_vars);
    auto noStates = (int)sylvan_satcount(states, s_vars);

    // mapper will store for each [decoded] state in the symbolic game, the parity game node id
    std::map<int, int> mapper;

    // Start constructing the parity game for given number of states (auto-grows)
    auto pargame = std::make_unique<pg::Game>(noStates);

    // index of next parity game vertex
    auto vidx = 0;

    // First initialize a parity game vertex for every state in the symbolic game
    {
        uint8_t state_arr[this->statebits];
        auto lf = mtbdd_enum_all_first(states, s_vars, state_arr, nullptr);
        while (lf != mtbdd_false) {
            // decode state
            auto state = 0;
            for (auto i=0; i<this->statebits; i++) {
                state <<= 1;
                if (state_arr[i]) state |= 1;
            }

            // create vertex controlled by Odd
            pargame->init_vertex(vidx, 0, 1, std::to_string(state));
            mapper[state] = vidx++;

            lf = mtbdd_enum_all_next(states, s_vars, state_arr, nullptr);
        }

        assert(vidx == noStates);
    }

    std::map<MTBDD, int> uap_to_vertex; // map after-uap to vertex
    std::map<MTBDD, int> cap_to_vertex; // map after-cap (priostate) to vertex

    uint8_t state_arr[this->statebits];
    auto lf = mtbdd_enum_all_first(this->trans, s_vars, state_arr, nullptr);
    while (lf != mtbdd_false) {
        // decode state
        auto state_i = 0;
        for (auto i=0; i<this->statebits; i++) {
            state_i <<= 1;
            if (state_arr[i]) state_i |= 1;
        }

        // translate state_i to state
        auto state_v = mapper.at(state_i);

        // find all "successors" of this state after the environment plays (uap)
        std::set<MTBDD> after_uap;
        {
            uint8_t uap_arr[this->uap_count];
            auto lf2 = mtbdd_enum_first(lf, uap_vars, uap_arr, nullptr);
            while (lf2 != mtbdd_false) {
                after_uap.insert(lf2);
                lf2 = mtbdd_enum_next(lf, uap_vars, uap_arr, nullptr);
            }
        }

        // start adding edges!
        pargame->e_start(state_v);

        for (auto uap_bdd : after_uap) {
            // check if we have seen this symbolic state before
            auto search = uap_to_vertex.find(uap_bdd);
            if (search != uap_to_vertex.end()) {
                // already exists
                // add edge from <state> to <search->second>
                pargame->e_add(state_v, search->second);
            } else {
                int uapv = vidx++;
                pargame->init_vertex(uapv, 0, 0, std::to_string(uap_bdd)); // controlled by Even
                uap_to_vertex[uap_bdd] = uapv;
                vertex_to_bdd[uapv] = uap_bdd;
                // add edge from <state> to <uapv>
                pargame->e_add(state_v, uapv);
            }
        }

        // we're done adding edges!
        pargame->e_finish();

        lf = mtbdd_enum_all_next(this->trans, s_vars, state_arr, nullptr);
    }

    // we now have all post-UAP in the uap_to_vertex map, so lets process them...

    for (auto & x : uap_to_vertex) {
        const MTBDD uap_bdd = x.first;
        const int uap_v = x.second;

        std::set<MTBDD> after_cap;
        {
            uint8_t cap_arr[this->cap_count];
            auto lf2 = mtbdd_enum_first(uap_bdd, cap_vars, cap_arr, nullptr);
            while (lf2 != mtbdd_false) {
                after_cap.insert(lf2);
                lf2 = mtbdd_enum_next(uap_bdd, cap_vars, cap_arr, nullptr);
            }
        }

        pargame->e_start(uap_v);

        for (auto cap_bdd : after_cap) {
            // only if not yet there
            auto search2 = cap_to_vertex.find(cap_bdd);
            if (search2 != cap_to_vertex.end()) {
                pargame->e_add(uap_v, search2->second);
            } else {
                int cap_v = vidx++;

                uint8_t pns_arr[this->priobits+this->statebits];
                mtbdd_enum_all_first(cap_bdd, pns_vars, pns_arr, nullptr);

                // decode priority of priostate
                int pr = 0;
                for (int i=0; i<this->priobits; i++) {
                    pr <<= 1;
                    if (pns_arr[i]) pr |= 1;
                }

                pargame->init_vertex(cap_v, pr, 1, std::to_string(cap_bdd)); // controlled by any
                cap_to_vertex[cap_bdd] = cap_v;
                vertex_to_bdd[cap_v] = cap_bdd;
                pargame->e_add(uap_v, cap_v);
            }
        }

        pargame->e_finish();
    }

    // we now have all post-CAP in the cap_to_vertex map, so lets process them...

    for (auto & x : cap_to_vertex) {
        const MTBDD cap_bdd = x.first;
        const int cap_v = x.second;

        uint8_t pns_arr[this->priobits+this->statebits];
        mtbdd_enum_all_first(cap_bdd, pns_vars, pns_arr, nullptr);

        // decode target
        int to = 0;
        for (int i=0; i<this->statebits; i++) {
            to <<= 1;
            if (pns_arr[this->priobits+i]) to |= 1;
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


/**
 * Restrict the transition relation to the given complete strategy.
 * Assumes <str> is defined on state + uap
 */
MTBDD select_str(MTBDD trans, MTBDD str)
{
    if (trans == str) {
        return mtbdd_true;
    } else if (mtbdd_isleaf(trans)) {
        return mtbdd_false;
    } else {
        // state or uncontrollable ap
        auto low = select_str(mtbdd_getlow(trans), str);
        mtbdd_refs_push(low);
        auto high = select_str(mtbdd_gethigh(trans), str);
        mtbdd_refs_push(high);
        MTBDD res = mtbdd_makenode(mtbdd_getvar(trans), low, high);
        mtbdd_refs_pop(2);
        return res;
    }
}


/**
 * Apply a strategy to the transition relation resulting in the strategy BDD
 * Returns mtbdd_invalid if something is wrong.
 */
typedef const std::map<MTBDD,MTBDD> strmap;
TASK_3(MTBDD, apply_str, MTBDD, trans, strmap*, str, int, first_cap)
{
    if (mtbdd_isleaf(trans)) {
        return mtbdd_false; // no unsanctioned leaves!
    } else {
        auto var = mtbdd_getvar(trans);
        if (var < (unsigned)first_cap) {
            // state or uncontrollable ap
            auto low = CALL(apply_str, mtbdd_getlow(trans), str, first_cap);
            mtbdd_refs_push(low);
            auto high = CALL(apply_str, mtbdd_gethigh(trans), str, first_cap);
            mtbdd_refs_push(high);
            if (low == mtbdd_invalid || high == mtbdd_invalid) {
                mtbdd_refs_pop(2);
                return mtbdd_invalid;
            } else {
                auto res = mtbdd_makenode(var, low, high);
                mtbdd_refs_pop(2);
                return res;
            }
        } else {
            auto it = str->find(trans);
                // FIXME HERE: should return the CAP resulting in it, not the thing!!!
            if (it != str->end()) {
                return select_str(it->first, it->second);
            } else {
                return mtbdd_false; // remove vertices without strategy
            }
        }
    }
}


bool SymGame::applyStrategy(const std::map<MTBDD, MTBDD>& str)
{
    auto strategy = RUN(apply_str, this->trans, &str, mtbdd_set_first(cap_vars));
    if (strategy == mtbdd_invalid) {
        std::cout << "Unable to compute strategy" << std::endl;
        return false;
    }

    this->strategies = strategy;
    return true;
}


std::unique_ptr<pg::Game> SymGame::strategyToPG()
{
    MTBDD states = sylvan_project(this->strategies, s_vars);
    auto noStates = (long)sylvan_satcount(states, s_vars);

    auto pargame = std::make_unique<pg::Game>(noStates);

    std::map<int, int> mapper;
    {
        auto idx = 0;

        uint8_t state_arr[this->statebits];
        auto lf = mtbdd_enum_all_first(states, s_vars, state_arr, nullptr);
        while (lf != mtbdd_false) {
            // decode state
            int state = 0;
            for (int i=0; i<this->statebits; i++) {
                state <<= 1;
                if (state_arr[i]) state |= 1;
            }

            pargame->init_vertex(idx, 0, 1, std::to_string(state)); // controlled by the Odd
            mapper[state] = idx++;

            lf = mtbdd_enum_all_next(states, s_vars, state_arr, nullptr);
        }
    }

    // DO A THING
    auto full = sylvan_and(this->trans, this->strategies);
    // We only care about priority 0 priostates
    while (!mtbdd_isleaf(full) && mtbdd_getvar(full) < (unsigned) this->priobits) {
        full = mtbdd_getlow(full);
    }

    auto vidx = (int)noStates;

    std::vector<int> succ_state;  // for current state, the successors

    uint8_t state_arr[this->statebits];
    MTBDD lf = mtbdd_enum_all_first(full, s_vars, state_arr, nullptr);
    while (lf != mtbdd_false) {
        // decode state
        auto _s = 0;
        for (auto i=0; i<this->statebits; i++) {
            _s <<= 1;
            if (state_arr[i]) _s |= 1;
        }

        int state = mapper.at(_s);

        // std::cout << "strategies for state " << state << ": " << lf << std::endl;

        auto pns = sylvan_project(lf, pns_vars);
        mtbdd_refs_pushptr(&pns);

        uint8_t pns_arr[this->priobits+this->statebits];
        auto lf2 = mtbdd_enum_all_first(pns, pns_vars, pns_arr, nullptr);
        while (lf2 != mtbdd_false) {
            assert(lf2 == mtbdd_true);

            // decode priostate
            int pr = 0;
            for (int i=0; i<this->priobits; i++) {
                pr <<= 1;
                if (pns_arr[i]) pr |= 1;
            }
            int to = 0;
            for (int i=0; i<this->statebits; i++) {
                to <<= 1;
                if (pns_arr[this->priobits+i]) to |= 1;
            }
            to = mapper.at(to);

            pargame->init_vertex(vidx, pr, 1);
            pargame->e_start(vidx);
            pargame->e_add(vidx, to);
            pargame->e_finish();
            succ_state.push_back(vidx);
            vidx++;

            // std::cout << "to (" << pr << ") " << to << std::endl;

            lf2 = mtbdd_enum_all_next(pns, pns_vars, pns_arr, nullptr);
        }

        pargame->e_start(state);
        for (int v : succ_state) pargame->e_add(state, v);
        pargame->e_finish();
        succ_state.clear();

        mtbdd_refs_popptr(1);

        lf = mtbdd_enum_all_next(full, s_vars, state_arr, nullptr);
    }

    mtbdd_refs_popptr(2);

    // tell Oink we're done adding stuff, resize game to final size
    pargame->v_resize(vidx);
    return pargame;
}


bool SymGame::solve(bool verbose)
{
    // prepare Odd (all odd-priority states)

    auto odd = sylvan_ithvar(this->priobits-1); // deepest bit is the parity
    mtbdd_refs_pushptr(&odd);

    // prepare renamers s_to_ns and ns_to_s

    auto pns_to_ps = mtbdd_map_empty();
    auto ps_to_pns = mtbdd_map_empty();
    mtbdd_refs_pushptr(&pns_to_ps);
    mtbdd_refs_pushptr(&ps_to_pns);
    {
        auto _s = ps_vars;
        auto _n = pns_vars;
        while (!mtbdd_set_isempty(_s)) {
            auto sv = mtbdd_set_first(_s);
            auto nv = mtbdd_set_first(_n);
            _s = mtbdd_set_next(_s);
            _n = mtbdd_set_next(_n);
            pns_to_ps = mtbdd_map_add(pns_to_ps, nv, sylvan_ithvar(sv));
            ps_to_pns = mtbdd_map_add(ps_to_pns, sv, sylvan_ithvar(nv));
        }
        assert(mtbdd_set_isempty(_n));
    }

    // prepare for every priority, priostates := the states of that priority

    MTBDD priostates[this->maxprio+1];
    for (int i=0; i<=this->maxprio; i++) {
        priostates[i] = encode_prio(i, this->priobits);
        mtbdd_refs_pushptr(&priostates[i]);
    }

    // prepare for every priority, lowereq := the states of <= priority

    MTBDD lowereq[this->maxprio+1];
    for (int i=0; i<=this->maxprio; i++) {
        lowereq[i] = mtbdd_false;
        mtbdd_refs_pushptr(&lowereq[i]);
        for (int j=0; j<=i; j++) {
            lowereq[i] = sylvan_or(priostates[j], lowereq[i]);
        }
    }

    // using the freezing FPI algorithm
    // prepare the distractions set
    MTBDD distractions = mtbdd_false; // initially, no distractions
    mtbdd_refs_pushptr(&distractions);

    // prepare the strategies_bdd set
    MTBDD strategies_bdd = mtbdd_false; // will have strategies_bdd ([prio]state -> priostate)
    mtbdd_refs_pushptr(&strategies_bdd);

    // prepare for every priority, the set of frozen states
    MTBDD freeze[this->maxprio+1];
    for (int i=0; i<=this->maxprio; i++) {
        freeze[i] = mtbdd_false;
        mtbdd_refs_pushptr(&freeze[i]);
    }

    int pr = 0;
    int iterations = 0;
    while (pr <= this->maxprio) {
        ++iterations;
        if (verbose) {
            // too verbose... comment out for now
            // std::cerr << "priority " << pr << std::endl;
        }

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
            onestepeven = sylvan_compose(onestepeven, ps_to_pns);
        }

        // then take product with transition
        // remember the strat: state -> uap -> cap
        MTBDD strat = onestepeven = sylvan_and_exists(this->trans, onestepeven, pns_vars);
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
            for (int i=pr+1; i<=this->maxprio; i++) {
                // remove each higher freeze set from <lwr>
                lwr = sylvan_and(lwr, sylvan_not(freeze[i]));
            }

            // HANDLING THE STRATEGY
            //   strat := good_state -> uap -> cap
            // (only good/useful strategies_bdd) this way, we only store useful strategies_bdd
            // first restrict strat to <lwr> (<=pr and not frozen higher)
            strat = sylvan_and(strat, lwr);
            // next restrict strat to not frozen lower (so we preserve the correct strategy)
            for (int i=0; i<=pr; i++) {
                strat = sylvan_and(strat, sylvan_not(freeze[i]));
            }
            // now strat is only for unfrozen vertices in the <=pr game

            // Extract strategy (for controller) -- str_vars is all except priorities
            // strat = sylvan_project(strat, str_vars);
            MTBDD states_in_strat = sylvan_project(strat, ps_vars);
            mtbdd_refs_pushptr(&states_in_strat);
            strategies_bdd = sylvan_ite(states_in_strat, strat, strategies_bdd); // big updater...
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
            for (int i=0; i<=this->maxprio; i++) {
                strat = sylvan_and(strat, sylvan_not(freeze[i]));
            }
            // Now strat is only for unfrozen vertices in the <=pr game
            // Update <strategies_bdd> with <strat> but only for states in <strat>
            MTBDD states_in_strat = sylvan_project(strat, ps_vars);
            mtbdd_refs_pushptr(&states_in_strat);
            strategies_bdd = sylvan_ite(states_in_strat, strat, strategies_bdd); // big updater...
            mtbdd_refs_popptr(1); // pop states_in_strat

            pr++;
        }

        mtbdd_refs_popptr(3); // pop newd, strat, onestepeven
    }

    if (verbose) {
        std::cerr << "symbolic solver required " << iterations << " iterations." << std::endl;
    }

    // we now know if the initial state is distracting or not
    // We force initial state to be state 0 in construction
    auto initial = encode_priostate(0, 0, s_vars, p_vars);
    mtbdd_refs_pushptr(&initial);

    if (sylvan_and(initial, distractions) != sylvan_false) {
        mtbdd_refs_popptr(6+3*this->maxprio); // free it up

        return false;
    }

    // We only care about priority 0 priostates
    while (!mtbdd_isleaf(strategies_bdd) && mtbdd_getvar(strategies_bdd) < (unsigned) this->priobits) {
        strategies_bdd = mtbdd_getlow(strategies_bdd);
    }

    mtbdd_refs_popptr(6+3+3*this->maxprio); // WHATEVER

    this->strategies = strategies_bdd;
    return true;
}


void SymGame::postprocess(bool verbose)
{
    // Select "lowest" strategy [heuristic]
    // TODO it would be nicer if we could do bisimulation without this heuristic
    strategies = RUN(clarify, strategies, cap_vars);

    // Now remove all unreachable states according to the strategy  (slightly smaller controller)
    {
        auto ns_to_s = mtbdd_map_empty();
        mtbdd_refs_pushptr(&ns_to_s);

        {
            auto _s = s_vars;
            auto _n = ns_vars;
            while (!mtbdd_set_isempty(_s)) {
                auto sv = mtbdd_set_first(_s);
                auto nv = mtbdd_set_first(_n);
                _s = mtbdd_set_next(_s);
                _n = mtbdd_set_next(_n);
                ns_to_s = mtbdd_map_add(ns_to_s, nv, sylvan_ithvar(sv));
            }
        }

        auto T = mtbdd_false;
        mtbdd_refs_pushptr(&T);

        {
            // compute the strategy transition relation
            // s > ns
            auto vars = np_vars;
            mtbdd_refs_pushptr(&vars);
            vars = mtbdd_set_addall(vars, cap_vars);
            vars = mtbdd_set_addall(vars, uap_vars);
            T = mtbdd_and_exists(strategies, this->trans, vars);
            mtbdd_refs_popptr(1);
        }

        auto visited = encode_state(0, s_vars);
        mtbdd_refs_pushptr(&visited);

        auto old = mtbdd_false;
        mtbdd_refs_pushptr(&old);

        // Just quick and dirty symbolic reachability ... this is not super efficient but it's OK
        while (old != visited) {
            old = visited;
            auto next = mtbdd_and_exists(visited, T, s_vars); // cross product
            mtbdd_refs_push(next);
            next = sylvan_compose(next, ns_to_s);   // and rename
            visited = sylvan_or(visited, next);
            mtbdd_refs_pop(1);
        }

        assert(visited == mtbdd_true || mtbdd_getvar(visited) >= (unsigned)this->priobits); // should not have priorities

        // count the number of states
        if (verbose) {
            auto noStates = (long)sylvan_satcount(visited, s_vars);
            std::cerr << "after reachability: " << noStates << " states." << std::endl;
        }

        strategies = sylvan_and(strategies, visited);
        trans = sylvan_and(trans, strategies);
        mtbdd_refs_popptr(4);
    }
}


[[maybe_unused]] void SymGame::print_vars() const
{
    {
        std::cerr << "s vars:";
        MTBDD _vars = s_vars;
        while (_vars != mtbdd_set_empty()) {
            std::cerr << " " << mtbdd_set_first(_vars);
            _vars = mtbdd_set_next(_vars);
        }
        std::cerr << std::endl;
    }

    {
        std::cerr << "u vars:";
        MTBDD _vars = uap_vars;
        while (_vars != mtbdd_set_empty()) {
            std::cerr << " " << mtbdd_set_first(_vars);
            _vars = mtbdd_set_next(_vars);
        }
        std::cerr << std::endl;
    }

    {
        std::cerr << "c vars:";
        MTBDD _vars = cap_vars;
        while (_vars != mtbdd_set_empty()) {
            std::cerr << " " << mtbdd_set_first(_vars);
            _vars = mtbdd_set_next(_vars);
        }
        std::cerr << std::endl;
    }

    {
        std::cerr << "p vars:";
        MTBDD _vars = np_vars;
        while (_vars != mtbdd_set_empty()) {
            std::cerr << " " << mtbdd_set_first(_vars);
            _vars = mtbdd_set_next(_vars);
        }
        std::cerr << std::endl;
    }

    {
        std::cerr << "ns vars:";
        MTBDD _vars = ns_vars;
        while (_vars != mtbdd_set_empty()) {
            std::cerr << " " << mtbdd_set_first(_vars);
            _vars = mtbdd_set_next(_vars);
        }
        std::cerr << std::endl;
    }
}


void SymGame::print_kiss(bool only_strategy)
{
    MTBDD _trans = this->trans;
    mtbdd_protect(&_trans);
    if (only_strategy) {
        _trans = sylvan_and(_trans, this->strategies);
    }
    MTBDD vars = mtbdd_set_empty();
    mtbdd_protect(&vars);
    // vars = mtbdd_set_addall(vars, this->s_vars);
    vars = mtbdd_set_addall(vars, this->uap_vars);
    vars = mtbdd_set_addall(vars, this->cap_vars);
    vars = mtbdd_set_addall(vars, this->pns_vars);
    std::stringstream ss;
    std::vector<std::string> lines;
    uint8_t s_arr[statebits];
    uint8_t arr[mtbdd_set_count(vars)+1];
    int states = 0;
    MTBDD slf = mtbdd_enum_all_first(_trans, s_vars, s_arr, nullptr);
    while (slf != mtbdd_false) {
        states++;
        // decode state
        int state = 0;
        for (int i=0; i<statebits; i++) {
            state <<= 1;
            if (s_arr[i]) state |= 1;
        }
        // find all things
        MTBDD lf = mtbdd_enum_first(slf, vars, arr, nullptr);
        while (lf != mtbdd_false) {
            int idx=0;
            std::string uap;
            // decode uap
            for (int i=0; i<this->uap_count; i++) {
                if (arr[idx] == 0) uap += "0";
                else if (arr[idx] == 1) uap += "1";
                else if (arr[idx] == 2) uap += "-";
                idx++;
            }
            // decode cap
            std::string cap;
            for (int i=0; i<this->cap_count; i++) {
                if (arr[idx] == 0) cap += "0";
                else if (arr[idx] == 1) cap += "1";
                else if (arr[idx] == 2) cap += "-";
                idx++;
            }
            // decode prio
            int prio = 0;
            for (int i=0; i<this->priobits; i++) {
                prio <<= 1;
                if (arr[idx++]) prio |= 1;
            }
            // decode next state
            int next_state = 0;
            for (int i=0; i<this->statebits; i++) {
                next_state <<= 1;
                if (arr[idx] == 1) next_state |= 1;
                else if (arr[idx] == 2) {
                    std::cerr << "ERROR: did not expect multiple next states!" << std::endl;
                }
                idx++;
            }
            assert(idx == mtbdd_set_count(vars));
            ss.str("");
            ss << uap << " s" << state << " s" << next_state << " " << cap;
            lines.push_back(ss.str());
            lf = mtbdd_enum_next(slf, vars, arr, nullptr);
        }
        // next state
        slf = mtbdd_enum_all_next(_trans, s_vars, s_arr, nullptr);
    }
    mtbdd_unprotect(&vars);
    mtbdd_unprotect(&_trans);

    std::cout << ".i " << uap_count << std::endl;
    std::cout << ".o " << cap_count << std::endl;
    std::cout << ".p " << lines.size() << std::endl;
    std::cout << ".s " << states << std::endl;
    std::cout << ".r 0" << std::endl;
    for (auto& line : lines) std::cout << line << std::endl;
}


[[maybe_unused]] void SymGame::print_trans(bool only_strategy) const
{
    MTBDD vars = mtbdd_set_empty();
    mtbdd_protect(&vars);
    vars = mtbdd_set_addall(vars, this->s_vars);
    vars = mtbdd_set_addall(vars, this->uap_vars);
    vars = mtbdd_set_addall(vars, this->cap_vars);
    vars = mtbdd_set_addall(vars, this->pns_vars);
    uint8_t arr[mtbdd_set_count(vars)+1];
    MTBDD _trans = this->trans;
    mtbdd_protect(&_trans);
    if (only_strategy) {
        _trans = sylvan_and(_trans, this->strategies);
    }
    MTBDD lf = mtbdd_enum_all_first(_trans, vars, arr, nullptr);
    while (lf != mtbdd_false) {
        int idx=0;
        // decode state
        int state = 0;
        for (int i=0; i<this->statebits; i++) {
            state <<= 1;
            if (arr[idx++]) state |= 1;
        }
        // decode uap
        int uap = 0;
        for (int i=0; i<this->uap_count; i++) {
            uap <<= 1;
            if (arr[idx++]) uap |= 1;
        }
        // decode cap
        int cap = 0;
        for (int i=0; i<this->cap_count; i++) {
            cap <<= 1;
            if (arr[idx++]) cap |= 1;
        }
        // decode prio
        int prio = 0;
        for (int i=0; i<this->priobits; i++) {
            prio <<= 1;
            if (arr[idx++]) prio |= 1;
        }
        // decode next state
        int next_state = 0;
        for (int i=0; i<this->statebits; i++) {
            next_state <<= 1;
            if (arr[idx++]) next_state |= 1;
        }
        assert(idx == mtbdd_set_count(vars));
        std::cerr << "from " << state << " uap " << uap << ": cap " << cap << " --> (" << prio << ") " << next_state << std::endl;
        lf = mtbdd_enum_all_next(_trans, vars, arr, nullptr);
    }
    mtbdd_unprotect(&vars);
    mtbdd_unprotect(&_trans);
}


[[maybe_unused]] void SymGame::print_strategies() const
{
    // strategy: s > uap > cap
    MTBDD vars = mtbdd_set_empty();
    mtbdd_protect(&vars);
    vars = mtbdd_set_addall(vars, this->s_vars);
    vars = mtbdd_set_addall(vars, this->uap_vars);
    vars = mtbdd_set_addall(vars, this->cap_vars);
    uint8_t arr[mtbdd_set_count(vars)+1];
    MTBDD lf = mtbdd_enum_all_first(this->strategies, vars, arr, nullptr);
    while (lf != mtbdd_false) {
        int idx=0;
        // decode state
        int state = 0;
        for (int i=0; i<this->statebits; i++) {
            state <<= 1;
            if (arr[idx++]) state |= 1;
        }
        // decode uap
        int uap = 0;
        for (int i=0; i<this->uap_count; i++) {
            uap <<= 1;
            if (arr[idx++]) uap |= 1;
        }
        // decode cap
        int cap = 0;
        for (int i=0; i<this->cap_count; i++) {
            cap <<= 1;
            if (arr[idx++]) cap |= 1;
        }
        assert(idx == mtbdd_set_count(vars));
        std::cerr << "from " << state << " uap " << uap << ": " << cap << std::endl;
        lf = mtbdd_enum_all_next(this->strategies, vars, arr, nullptr);
    }
    mtbdd_unprotect(&vars);
}


sylvan::MTBDD SymGame::evalLabel(BTree* label, HoaData* data, uint32_t* variables) {
    return RUN(evalLabel, label, data, variables);
}