/**
 * Copyright Tom van Dijk
 */

#ifndef KNOR_PERMISSIVE_HPP
#define KNOR_PERMISSIVE_HPP

#include <vector>

#include <oink/game.hpp>

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
    };

    /**
     * The positional winning strategy that Oink computed, as a multi-strategy:
     * the successor of each vertex owned and won by Even (empty for others).
     */
    static std::vector<std::vector<int>> positional(const pg::Game& game);

    /**
     * Solve the game by fixpoint iteration with freezing (FPI), as the
     * symbolic solver does, and return its multi-strategy: for each vertex
     * owned and won by Even, every successor that is good for Even when the
     * vertex is decided (empty for other vertices). Any choice among them is
     * the strategy of FPI for some order of the edges, hence winning. Sets
     * winner to the winner (0 or 1) of each vertex.
     */
    static std::vector<std::vector<int>> fpim(const pg::Game& game, std::vector<int>& winner);

    /**
     * The same with fixpoint iteration with justifications (FPJ): the moves of
     * a vertex are those good for its owner when it is justified, and a
     * justified vertex is only invalidated when its last move is.
     */
    static std::vector<std::vector<int>> fpjm(const pg::Game& game, std::vector<int>& winner);

    /**
     * Greedily extend a winning multi-strategy (as given by positional(),
     * fpim() or fpjm()): visit the vertices reachable from the initial vertex
     * breadth-first, and allow each further successor of a vertex of Even if
     * that keeps the multi-strategy winning, i.e., if the new edge closes no
     * cycle whose maximal priority is odd. Returns the allowed successors of
     * each vertex owned and won by Even (empty for others).
     */
    static std::vector<std::vector<int>> extend(const pg::Game& game, int initial,
                                                std::vector<std::vector<int>> allowed, Stats& stats);

    /**
     * Check independently that the multi-strategy is winning from the initial
     * vertex: no cycle reachable under it has an odd maximal priority.
     */
    static bool verify(const pg::Game& game, int initial, const std::vector<std::vector<int>>& allowed);
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
