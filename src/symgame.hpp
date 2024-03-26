/**
 * Copyright Tom van Dijk
 */

#include <sylvan.h>
#include <oink/oink.hpp>
#include <memory>

extern "C" {
    #include "simplehoa.h"
}

#ifndef KNOR_SYMGAME_HPP
#define KNOR_SYMGAME_HPP

class SymGame final {
public:
    int maxprio;

    sylvan::MTBDD p_vars;
    sylvan::MTBDD s_vars;
    sylvan::MTBDD uap_vars;
    sylvan::MTBDD cap_vars;
    sylvan::MTBDD np_vars;
    sylvan::MTBDD ns_vars;
    sylvan::MTBDD ps_vars;
    sylvan::MTBDD pns_vars;
    sylvan::MTBDD uns_vars;

    int statebits;
    int priobits;
    int cap_count;
    int uap_count;

    sylvan::MTBDD trans;       // transition relation of symbolic game:        state -> uap -> cap -> priority -> next_state
    sylvan::MTBDD strategies;  // contains the solution: the strategies:       good_state -> uap -> cap

    SymGame(int statebits, int priobits, int uap_count, int cap_count, int maxprio);
    SymGame(const SymGame&) = delete;
    ~SymGame();

    static std::unique_ptr<SymGame> constructSymGame(HoaData* data, bool isMaxParity, bool controllerIsOdd);

    /**
     * Encode a state as a BDD, using statebits 0..<statebits>, offsetted by <offset>+<priobits>
     * High-significant bits come before low-significant bits in the BDD
     */
    static sylvan::MTBDD encode_state(uint32_t state, sylvan::MTBDD state_vars);

    /**
     * Encode priority i.e. all states via priority <priority>
     */
    static sylvan::MTBDD encode_prio(int priority, int priobits);

    /**
     * Encode a priostate as a BDD, using priobits before statebits
     * High-significant bits come before low-significant bits in the BDD
     */
    static sylvan::MTBDD encode_priostate(uint32_t state, uint32_t priority, sylvan::MTBDD statevars, sylvan::MTBDD priovars);

    /**
     * Translate symbolic PG to explicit game in Oink, that can then be solved.
     */
    std::unique_ptr<pg::Game> toExplicit(std::map<int, sylvan::MTBDD>& vertex_to_bdd);

    /**
     * Convert the strategy of a realizable parity game to an explicit parity game that
     * should be won by player Even.
     */
    std::unique_ptr<pg::Game> strategyToPG();

    /**
     * Apply a strategy (computed via Oink)
     */
    bool applyStrategy(const std::map<sylvan::MTBDD, sylvan::MTBDD>&);

    /**
     * Solve the symbolic parity game, return true if won
     */
    bool solve(bool verbose);

    /**
     * After solving the game, compute BDDs for output and state
     */
    void postprocess(bool verbose);

    [[maybe_unused]] void print_vars() const;
    [[maybe_unused]] void print_kiss(bool only_strategy=false);
    [[maybe_unused]] void print_trans(bool only_strategy=false) const;
    [[maybe_unused]] void print_strategies() const;
};

#endif // KNOR_SYMGAME_HPP
