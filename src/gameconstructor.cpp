/**************************************************************************
 * Copyright Tom van Dijk
 *************************************************************************/

#include <cassert> // for assert
#include <cstring> // for strcmp
#include <iostream>
#include <set>

#include <oink/game.hpp>
#include <oink/solvers.hpp>
#include <symgame.hpp> // for encode_priostate
#include <bddtools.hpp>
#include <gameconstructor.hpp>

extern "C" {
#include <sylvan.h>
#include "simplehoa.h"
}

using namespace sylvan;


/**
 * Given a label and a valuation of some of the atomic propositions,
 * we determine whether the label is true (1), false (-1), or its
 * value is unknown (0). The valuation is expected as an unsigned
 * integer whose i-th bit is 1 iff the i-th AP in apIds is set to 1
 * Returns -2 if there is a problem.
 */
static int
evalLabelNaive(BTree* label, Alias* aliases, int numAliases, int numAPs, int* apIds, uint64_t value) {
    assert(label != nullptr);
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
 * This helper function ensures that the priority p is adjusted to
 * ensure we have a "max, even" parity game, as this is what Oink expects.
 * @param controllerIsOdd if the controller is the odd player
 * @param noPriorities how many priorities are in the game
 * @param maxPriority if the game is a max game
 */
int GameConstructor::adjustPriority(int p, bool maxPriority, bool controllerIsOdd, int noPriorities)
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
 * Construct and solve the game explicitly
 */
TASK_3(pg::Game*, constructGameNaive, HoaData*, data, bool, isMaxParity, bool, controllerIsOdd)
{
    // Set which APs are controllable in the bitset controllable
    pg::bitset controllable(data->noAPs);
    for (int i=0; i<data->noCntAPs; i++) {
        controllable[data->cntAPs[i]] = true;
    }

    // Count the number of controllable/uncontrollable APs
    const auto cap_count = controllable.count();
    const auto uap_count = data->noAPs - cap_count;
    const uint64_t numValuations = (1ULL << uap_count);

    int ucntAPs[uap_count];
    int uidx = 0;
    for (int i = 0; i < data->noAPs; i++) {
        if (!controllable[i]) ucntAPs[uidx++] = i;
    }

    // Now initialize a new parity game
    auto *game = new pg::Game(data->noStates * 10); // start with 10 the number of states
    // notice that the number of vertices automatically grows when needed anyway

    auto nextIndex = data->noStates; // index of the next vertex to make

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
                if (state->label != nullptr) label = state->label;
                else label = trans->label;
                assert(label != nullptr);
                // we add a vertex + edges if the transition is compatible with the
                // valuation we are currently considering
                int evald = evalLabelNaive(label, data->aliases, data->noAliases, (int)uap_count, ucntAPs, value);
                if (evald == -1) continue; // not compatible
                // there should be a priority at state or transition level
                if (state->accSig == nullptr) {
                    // there should be exactly one acceptance set!
                    assert(trans->noAccSig == 1);
                    // adjust priority
                    int priority = GameConstructor::adjustPriority(trans->accSig[0], isMaxParity, controllerIsOdd, data->noAccSets);

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

            auto vinter = nextIndex++;
            succ_state.push_back(vinter);
            game->init_vertex(vinter, 0, 0);
            game->e_start(vinter);
            for (auto to : succ_inter) game->e_add(vinter, to);
            game->e_finish();
            succ_inter.clear();
        }

        int priority;
        if (state->accSig != nullptr) {
            priority = GameConstructor::adjustPriority(state->accSig[0], isMaxParity, controllerIsOdd, data->noAccSets);
        } else {
            priority = 0;
        }

        game->init_vertex(state->id, priority, 1, state->name ? state->name : std::to_string(state->id));
        game->e_start(state->id);
        for (auto to : succ_state) game->e_add(state->id, to);
        game->e_finish();
        succ_state.clear();
    }

    // tell Oink we're done adding stuff, resize game to final size
    game->v_resize(nextIndex);

    return game;
}


/**
 * Given some intermediary MTBDD root, collect all the MTBDD leafs.
 * It decodes the leaf <priority, state> and re-encodes the leaf as a BDD such that
 * the first <priobits> BDD variables encode the priority and
 * the next <statebits> BDD variables encode the target state
 */
MTBDD collect_targets(MTBDD trans, std::set<uint64_t> &res, MTBDD statevars, MTBDD priovars)
{
    if (mtbdd_isleaf(trans)) {
        uint64_t leaf = mtbdd_getint64(trans);
        res.insert(leaf);

        auto priority = (uint32_t)(leaf>>32);
        auto state = (uint32_t)(leaf & 0xffffffff);

        return BDDTools::encode_priostate(state, priority, statevars, priovars);
    } else {
        auto left = mtbdd_false;
        auto right = mtbdd_false;
        mtbdd_refs_pushptr(&left);
        mtbdd_refs_pushptr(&right);

        left = collect_targets(mtbdd_getlow(trans), res, statevars, priovars);
        right = collect_targets(mtbdd_gethigh(trans), res, statevars, priovars);
        auto result = sylvan_or(left, right);

        mtbdd_refs_popptr(2);
        return result;
    }
}


/**
 * Construct and solve the game explicitly
 */
TASK_3(pg::Game*, constructGame, HoaData *, data, bool, isMaxParity, bool, controllerIsOdd)
{
    // Set which APs are controllable in the bitset controllable
    pg::bitset controllable(data->noAPs);
    for (int i=0; i<data->noCntAPs; i++) {
        controllable[data->cntAPs[i]] = true;
    }

    // Count the number of controllable/uncontrollable APs
    const auto cap_count = controllable.count();
    const auto uap_count = data->noAPs - cap_count;

    // Initialize the BDD variable indices
    // Variable order uncontrollable < controllable
    uint32_t variables[data->noAPs];
    uint32_t uidx = 0, oidx = data->noAPs - controllable.count();
    for (int i=0; i<data->noAPs; i++) {
        if (controllable[i]) variables[i] = oidx++;
        else variables[i] = uidx++;
    }

    // Now initialize a new parity game
    auto game = new pg::Game(data->noStates * 10); // start with 10 the number of states
    // The number of vertices automatically grows when needed, but 10xstates is a good start

    auto nextIndex = data->noStates; // index of the next vertex to make

    // Prepare number of statebits and priobits
    int statebits = 1;
    while ((1ULL<<(statebits)) <= (uint64_t)data->noStates) statebits++;
    int evenMax = 2 + 2*((data->noAccSets+1)/2); // should be enough...
    int priobits = 1;
    while ((1ULL<<(priobits)) <= (unsigned)evenMax) priobits++;

    // Prepare the p_vars and s_vars (for building the BDD)
    auto p_vars = mtbdd_set_empty();
    auto s_vars = mtbdd_set_empty();
    mtbdd_refs_pushptr(&p_vars);
    mtbdd_refs_pushptr(&s_vars);
    auto var = 0;
    for (auto i=0; i<priobits; i++) p_vars = mtbdd_set_add(p_vars, var++);
    for (auto i=0; i<statebits; i++) s_vars = mtbdd_set_add(s_vars, var++);

    std::vector<int> succ_state;  // for current state, the successors
    std::vector<int> succ_inter;  // for current intermediate state, the successors

    auto trans_bdd = mtbdd_false;
    mtbdd_refs_pushptr(&trans_bdd);

    std::set<uint64_t> targets;

    std::map<MTBDD, int> inter_vertices;
    std::map<uint64_t, int> target_vertices;

    // these variables are used deeper in the for loops, however we can
    // push them to mtbdd refs here and avoid unnecessary pushing and popping
    auto lblbdd = mtbdd_false;
    auto leaf = mtbdd_false;
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
            if (state->label != nullptr) label = state->label;
            else label = trans->label;
            assert(label != nullptr);
            // there should be a priority at state or transition level
            int priority = 0;
            if (state->accSig == nullptr) {
                // there should be exactly one acceptance set!
                assert(trans->noAccSig == 1);
                // adjust priority
                priority = GameConstructor::adjustPriority(trans->accSig[0], isMaxParity, controllerIsOdd, data->noAccSets);
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

        // TODO: this assumes uap_count is the first non-uap variable index
        auto inter_bdds = BDDTools::collectSubroots(trans_bdd, uap_count);
        for (MTBDD inter_bdd : inter_bdds) {
            auto targets_bdd = collect_targets(inter_bdd, targets, s_vars, p_vars);

            int vinter;
            auto it = inter_vertices.find(targets_bdd);
            if (it == inter_vertices.end()) {
                for (uint64_t lval : targets) {
                    int priority = (int)(lval >> 32);
                    int target = (int)(lval & 0xffffffff);

                    if (priority != 0) {
                        int vfin;
                        auto it2 = target_vertices.find(lval);
                        if (it2 == target_vertices.end()) {
                            vfin = nextIndex++;
                            game->init_vertex(vfin, priority, 0);
                            game->e_start(vfin);
                            game->e_add(vfin, target);
                            game->e_finish();
                            target_vertices.insert(std::make_pair(lval, vfin));
                        } else {
                            vfin = it2->second;
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
            priority = GameConstructor::adjustPriority(state->accSig[0], isMaxParity, controllerIsOdd, data->noAccSets);
        } else {
            priority = 0;
        }

        game->init_vertex(state->id, priority, 1, state->name ? state->name : std::to_string(state->id));
        game->e_start(state->id);
        for (int to : succ_state) game->e_add(state->id, to);
        game->e_finish();

        succ_state.clear();
        inter_vertices.clear();
    }

    // tell Oink we're done adding stuff, resize game to final size
    game->v_resize(nextIndex);

    mtbdd_refs_popptr(5);
    mtbdd_refs_pop(ref_counter);

    return game;
}


/**
 * Construct the symbolic game
 */
std::unique_ptr<SymGame> GameConstructor::constructSymGame(HoaData *data, bool isMaxParity, bool controllerIsOdd) {
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
                priority = GameConstructor::adjustPriority(trans->accSig[0], isMaxParity, controllerIsOdd, data->noAccSets);
            } else {
                auto id = trans->successors[0];
                assert(data->states[id].noAccSig == 1);
                priority = GameConstructor::adjustPriority(data->states[id].accSig[0], isMaxParity, controllerIsOdd, data->noAccSets);
            }
            if (priority > res->maxprio) res->maxprio = priority;
            // swap with initial if needed
            int succ = trans->successors[0];
            if (succ == 0) succ = vstart;
            else if (succ == vstart) succ = 0;
            // encode the label as a MTBDD
            lblbdd = RUN(evalLabel, label, data, variables);
            // encode priostate (leaf) and update transition relation
            leaf = BDDTools::encode_priostate(succ, priority, res->ns_vars, res->np_vars);
            // trans := lbl THEN leaf ELSE trans
            transbdd = mtbdd_ite(lblbdd, leaf, transbdd);
            // deref lblbdd and leaf
            lblbdd = leaf = mtbdd_false;
        }

        // encode source state and add to full transition relation
        int src = state->id;
        if (src == 0) src = vstart;
        else if (src == vstart) src = 0;
        statebdd = BDDTools::encode_state(src, res->s_vars);
        // update full trans with <statebdd> then <transbdd>
        res->trans = mtbdd_ite(statebdd, transbdd, res->trans);
        // deref statebdd and transbdd
        statebdd = transbdd = mtbdd_false;
    }

    mtbdd_refs_popptr(4);

    return res;
}


pg::Game* GameConstructor::constructGameNaive(HoaData *data, bool isMaxParity, bool controllerIsOdd) {
    return RUN(constructGameNaive, data, isMaxParity, controllerIsOdd);
}


pg::Game* GameConstructor::constructGame(HoaData *data, bool isMaxParity, bool controllerIsOdd) {
    return RUN(constructGame, data, isMaxParity, controllerIsOdd);
}
