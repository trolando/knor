/**
 * Copyright Tom van Dijk
 */

#include <sylvan.h>
#include <oink/oink.hpp>

extern "C" {
    #include "simplehoa.h"
}

#pragma once

class SymGame {
public:
    int maxprio;

    sylvan::MTBDD s_vars;
    sylvan::MTBDD uap_vars;
    sylvan::MTBDD cap_vars;
    sylvan::MTBDD p_vars;
    sylvan::MTBDD ns_vars;
    sylvan::MTBDD pns_vars;

    int statebits;
    int priobits;
    int cap_count;
    int uap_count;

    sylvan::MTBDD trans;       // transition relation of symbolic game:        state -> uap -> cap -> priority -> next_state
    sylvan::MTBDD strategies;  // contains the solution: the strategies:       good_state -> uap -> cap

    SymGame(int statebits, int priobits, int uap_count, int cap_count, int maxprio);
    virtual ~SymGame() ;

    /**
     * Translate symbolic PG to explicit game in Oink, that can then be solved.
     */
    pg::Game *toExplicit(std::map<int, sylvan::MTBDD>&);

    /**
     * Convert the strategy of a realizable parity game to an explicit parity game that
     * should be won by player Even.
     */
    pg::Game *strategyToPG();

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

    void print_trans();
    void print_strategies();
};


TASK_DECL_3(SymGame*, constructSymGame, HoaData* /*data*/, bool /*isMaxParity*/, bool /*controllerIsOdd*/)
TASK_DECL_2(bool, wrap_solve, SymGame* /*game*/, bool /*verbose*/)
VOID_TASK_DECL_2(wrap_pp, SymGame* /*game*/, bool /*verbose*/)

