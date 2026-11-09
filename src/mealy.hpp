/**
 * Copyright Tom van Dijk
 */

#ifndef KNOR_MEALY_HPP
#define KNOR_MEALY_HPP

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <symgame.hpp>

/**
 * Minimisation of the controller as a Mealy machine.
 *
 * The (multi-)strategy of a solved game, restricted to the states that are
 * reachable from the initial state, is a nondeterministic Mealy machine: in
 * every state, for every input, it allows some outputs, each with a
 * successor state. Because the strategy is winning, every deterministic
 * Mealy machine that only makes allowed moves is a correct controller, so
 * the parity condition plays no role here. minimize() looks for such a
 * machine with few states: after merging bisimilar states, it computes which
 * states are compatible (can be served by one controller state) and a set of
 * pairwise incompatible states, a lower bound. These anchor the first states
 * of the controller, so that each state can only be served by the anchors it
 * is compatible with, which keeps the SAT encoding small. It first tries as
 * many states as the lower bound (then the result is minimal), then searches
 * upwards with growing steps and bisects once a controller is found.
 */
class MealyMinimizer final {
public:
    /**
     * One search for a small controller within a set of moves (for --stats).
     */
    struct Search {
        std::string moves;      // "strategy", "extended strategy" or "all winning moves"
        std::string skipped;    // why the search did not run, empty if it did
        int machine = 0;        // states of the Mealy machine of the moves
        int quotient = 0;       // after merging bisimilar states
        int lowerBound = 0;     // pairwise incompatible states
        int found = 0;          // states of the controller found (the quotient if no smaller one), 0 if none
        bool proven = false;    // found is minimal for these moves
        struct Call {
            int states;         // the number of controller states tried
            int outcome;        // 1 suffice, 0 do not suffice, -1 undecided, -2 too large (not tried)
            long work;          // see Options
        };
        std::vector<Call> calls;
    };

    /**
     * With relax: another game of the controller (other state codes, more moves) for --dc.
     */
    enum class Kind {
        controller,             // a controller of another search (not the one returned)
        relaxed,                // with every move that keeps the controller correct, numbered by them
        cubeCodes,              // relaxed, with state codes as cubes
        adjacencyCodes,         // relaxed, with codes that differ in few bits for related states
        outputCodes,            // with output bits as state bits
        monitorCodes,           // the controller's own moves, with monitor flags as state bits
        relaxedMonitorCodes,    // relaxed, with monitor flags as state bits
        unclarified,            // without a minimised controller: the freedom of the strategy itself
    };
    struct Variant {
        std::unique_ptr<SymGame> game;
        std::string name;       // what it is, e.g. "relaxed, breadth-first"
        Kind kind;
        bool chosen = true;     // of the controller that minimize() returns, not of another one
    };

    struct Result {
        std::unique_ptr<SymGame> game;  // the minimised controller, or nullptr if none was found
        int states = 0;                 // states of the minimised controller
        std::vector<Variant> variants;  // with relax
        std::vector<Search> searches;   // the searches that were made
    };

    /**
     * The variants that relax makes; all by default (knor's ablation options
     * switch them off).
     */
    struct Variants {
        bool relaxed = true;        // with every move that keeps the controller correct (else its own moves)
        bool renumber = true;       // relaxed: number the states breadth-first and depth-first by these moves
        bool outputCodes = true;    // state codes with output bits
        bool monitorCodes = true;   // state codes with monitor flags
        bool cubeCodes = true;      // relaxed: state codes as cubes, so no code is unused
        bool adjacencyCodes = true; // relaxed: related states get codes that differ in few bits
    };

    /**
     * The SAT searches are limited by their work, which follows the time
     * closely (the propagations of the SAT solver, the clauses it inspects and
     * the clauses added; about 25 million per second), not by time, so that the
     * result does not depend on the machine or its load.
     */
    struct Options {
        long work = 2000000000;                 // for the search within the strategy, and the extended one
        bool verbose = false;
        bool relax = false;
        Variants variants;                      // with relax
        sylvan::MTBDD extended = sylvan::mtbdd_false;   // a winning multi-strategy that contains the strategy
        sylvan::MTBDD anyMoves = sylvan::mtbdd_false;   // all moves into the winning region
        long anyWork = 2000000000;              // for the search over anyMoves
        bool fewerFirst = true;                 // after an undecided SAT call, try fewer states first
        bool strongBound = true;                // lower bound from all unavoidable states, the initial one optional
        bool winningBound = true;               // over anyMoves, from the states that every winning controller visits
        int anyStates = 0;                      // search anyMoves only for controllers of at most this many states (0: all)
    };

    /**
     * Minimise the controller given by game.trans and game.strategies, using
     * at most options.work. The result is a game with the controller's
     * states as states, in which strategies and trans contain exactly the
     * controller's moves (with priority 0). With relax, the
     * variants are games of the same controller with other state codes,
     * most with every move that keeps the controller correct (possibly
     * several outputs and successors per state and input); only AIGEncoder
     * with setDontCares() may encode them. Then it searches the wider move
     * sets of the options that are not mtbdd_false, each with at most its
     * work, for a smaller controller that may make any of their moves: first
     * extended, then anyMoves (only if the controller has at most anyStates
     * states; such controllers are checked for losing cycles unless the moves
     * form a winning multi-strategy). More freedom weakens the lower bound,
     * so this gives the best of both: the fast search of the strategy, and
     * smaller controllers where the wider search finds them. With relax, a
     * wider controller that is not smaller is a variant too.
     */
    static Result minimize(SymGame& game, const Options& options);
};

#endif
