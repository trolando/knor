/**
 * Copyright Tom van Dijk
 */

#include <set>
#include <sylvan.h>
#include <sylvan_mtbdd.h>

#ifndef KNOR_BDDTOOLS_HPP
#define KNOR_BDDTOOLS_HPP

class BDDTools {
public:
    /**
     * Find all BDD subroots with first variable >= firstVar
     * @param root the BDD to find the subroots of
     * @param firstVar first variable of subroots
     * @param collector set to collect subroots in
     */
    static void collectSubroots(sylvan::MTBDD root, uint32_t firstVar, std::set<sylvan::MTBDD> &collector);

    static std::set<sylvan::MTBDD> collectSubroots(sylvan::MTBDD bdd, uint32_t firstVar) {
        std::set<sylvan::MTBDD> res;
        collectSubroots(bdd, firstVar, res);
        return res;
    }

    /**
     * Find all paths from the given root to the given subroot.
     * @param root
     * @param subroot
     * @return a BDD with all paths to the given subroot.
     */
    static sylvan::MTBDD pathsToSubroot(sylvan::MTBDD root, uint32_t firstVar, sylvan::MTBDD subroot);
};

#endif //KNOR_BDDTOOLS_HPP
