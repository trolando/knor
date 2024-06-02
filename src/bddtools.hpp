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

    /**
     * Encode a priostate as a BDD, using priobits before statebits
     * High-significant bits come before low-significant bits in the BDD
     */
    static sylvan::MTBDD encode_priostate(uint32_t state, uint32_t priority, sylvan::MTBDD statevars, sylvan::MTBDD priovars);

    /**
     * Encode a state as a BDD, using statebits 0..<statebits>, offsetted by <offset>+<priobits>
     * High-significant bits come before low-significant bits in the BDD
     */
    static sylvan::MTBDD encode_state(uint32_t state, sylvan::MTBDD state_vars);

    /**
     * Encode priority i.e. all states via priority <priority>
     */
    static sylvan::MTBDD encode_prio(int priority, int priobits);

};

#endif //KNOR_BDDTOOLS_HPP
