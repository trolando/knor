/**
 * Copyright Tom van Dijk
 */

#include "bddtools.hpp"
#include <vector>

using namespace sylvan;

void BDDTools::collectSubroots(MTBDD root, uint32_t firstVar, std::set<sylvan::MTBDD> &collector)
{
    std::vector<MTBDD> stack;
    std::set<MTBDD> seen;
    stack.push_back(root);
    seen.insert(root);

    while (!stack.empty()) {
        MTBDD current = stack.back();
        stack.pop_back();
        if (current != mtbdd_false) {
            if (mtbdd_isleaf(current)) {
                collector.insert(current);
            } else {
                if (mtbdd_getvar(current) < firstVar) {
                    auto high = mtbdd_gethigh(current);
                    auto low = mtbdd_getlow(current);
                    if (seen.emplace(high).second) stack.push_back(high);
                    if (seen.emplace(low).second) stack.push_back(low);
                } else {
                    collector.insert(current);
                }
            }
        }
    }
}

TASK_3(MTBDD, paths_to_subroot, MTBDD, bdd, uint32_t, firstvar, MTBDD, subroot)
{
    // TODO might need caching...
    if (bdd == subroot) {
        return mtbdd_true;
    } else if (mtbdd_isleaf(bdd) || mtbdd_getvar(bdd) >= firstvar) {
        return mtbdd_false;
    } else {
        auto var = mtbdd_getvar(bdd);
        mtbdd_refs_spawn(SPAWN(paths_to_subroot, mtbdd_gethigh(bdd), firstvar, subroot));
        MTBDD low = mtbdd_refs_push(CALL(paths_to_subroot, mtbdd_getlow(bdd), firstvar, subroot));
        MTBDD high = mtbdd_refs_push(mtbdd_refs_sync(SYNC(paths_to_subroot)));
        MTBDD res = mtbdd_makenode(var, low, high);
        mtbdd_refs_pop(2);
        return res;
    }
}

MTBDD BDDTools::pathsToSubroot(MTBDD root, uint32_t firstVar, MTBDD subroot)
{
    return RUN(paths_to_subroot, root, firstVar, subroot);
}

/**
 * Encode a state as a BDD, using statebits 0..<statebits>, offsetted by <offset>+<priobits>
 * High-significant bits come before low-significant bits in the BDD
 */
MTBDD BDDTools::encode_state(uint32_t state, MTBDD statevars)
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
MTBDD BDDTools::encode_prio(int priority, int priobits)
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
MTBDD BDDTools::encode_priostate(uint32_t state, uint32_t priority, MTBDD statevars, MTBDD priovars)
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
