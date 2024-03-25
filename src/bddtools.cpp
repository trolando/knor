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
