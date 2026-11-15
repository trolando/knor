/**
 * Copyright Tom van Dijk
 */

#ifndef KNOR_PERMISSIVE_HPP
#define KNOR_PERMISSIVE_HPP

#include <vector>

#include <oink/game.hpp>
#include <oink/oink.hpp>
#include <oink/solution.hpp>

/**
 * Multi-strategies for player Even (the controller) in a parity game solved
 * by Oink (maximal priority wins, even priorities are good for Even). A
 * multi-strategy allows every vertex owned and won by Even a set of
 * successors; it is winning if every play from the initial vertex in which
 * Even only takes allowed edges is won by Even.
 */
class PermissiveStrategy final {
public:
    struct Stats {
        long vertices = 0;      // reachable vertices of Even that were extended
        long candidates = 0;    // edges that were considered
        long added = 0;         // edges that were added
        long work = 0;          // (vertex, level) pairs visited by the cycle checks
        bool stopped = false;   // the budget ran out
    };

    /**
     * Prepare Oink to solve with fixpoint iteration (fpim), or with
     * justifications (fpjm), and keep the multi-strategy that it records: for
     * each vertex, every move that is good for its owner when the vertex is
     * decided (or justified). Any choice among them is the strategy of the
     * solver for some order of the edges, hence winning. Oink's preprocessing
     * (self-loops, winner-controlled winning cycles, games with one parity) and
     * solving by bottom SCCs would decide some vertices with one move each, so
     * they are off.
     */
    static void multiSolver(pg::Oink& engine, bool justifications);

    /**
     * The strategy of a solution as a multi-strategy: for each vertex owned and
     * won by Even, its moves (the one of a positional strategy, or all of a
     * multi-strategy); empty for other vertices.
     */
    static std::vector<std::vector<int>> moves(const pg::Game& game, const pg::Solution& solution);

    /**
     * Greedily extend a winning multi-strategy (as given by moves()): visit
     * the vertices reachable from the initial vertex
     * breadth-first, and allow each further successor of a vertex of Even if
     * that keeps the multi-strategy winning, i.e., if the new edge closes no
     * cycle whose maximal priority is odd. With a budget (not 0), no further
     * edges are considered once the cycle checks have visited that many
     * (vertex, level) pairs; the result is still winning. Returns the allowed
     * successors of each vertex owned and won by Even (empty for others).
     */
    static std::vector<std::vector<int>> extend(const pg::Game& game, const pg::Solution& solution, int initial,
                                                std::vector<std::vector<int>> allowed, long budget, Stats& stats);

    /**
     * Check independently that the multi-strategy is winning from the initial
     * vertex: no cycle reachable under it has an odd maximal priority.
     */
    static bool verify(const pg::Game& game, const pg::Solution& solution, int initial,
                       const std::vector<std::vector<int>>& allowed);
};

/**
 * Is there a cycle reachable from the initial vertex whose maximal priority is
 * odd, in the graph with the given successors and vertex priorities?
 */
bool hasOddMaxCycle(const std::vector<std::vector<int>>& successors, const std::vector<int>& priorities, int initial);

/**
 * Such a cycle, as its vertices in order (the last one has an edge to the
 * first), or an empty list if there is none.
 */
std::vector<int> oddMaxCycle(const std::vector<std::vector<int>>& successors, const std::vector<int>& priorities, int initial);

#endif
