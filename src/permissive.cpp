/**
 * Copyright Tom van Dijk
 */

#include <algorithm>
#include <cstdint>
#include <deque>
#include <utility>
#include <vector>

#include <permissive.hpp>

namespace {

// Compress priorities to levels 0..d-1, keeping their order; odd[l] tells
// whether the priorities of level l are odd.
void toLevels(const std::vector<int>& priority, std::vector<int>& level, std::vector<bool>& odd)
{
    std::vector<int> distinct(priority);
    std::sort(distinct.begin(), distinct.end());
    distinct.erase(std::unique(distinct.begin(), distinct.end()), distinct.end());
    odd.resize(distinct.size());
    for (size_t i = 0; i < distinct.size(); i++) odd[i] = distinct[i] % 2 != 0;
    level.resize(priority.size());
    for (size_t v = 0; v < priority.size(); v++) {
        level[v] = static_cast<int>(std::lower_bound(distinct.begin(), distinct.end(), priority[v]) - distinct.begin());
    }
}

// The game as restricted by a multi-strategy: vertices owned and won by Even
// have their allowed successors, all other vertices all their successors.
class Restricted final {
public:
    Restricted(const pg::Game& game, const std::vector<std::vector<int>>& allowed)
        : game(game), allowed(allowed)
    {
        std::vector<int> priority(game.vertexcount());
        for (int v = 0; v < game.vertexcount(); v++) priority[v] = game.priority(v);
        toLevels(priority, level, odd);
    }

    int vertices() const { return static_cast<int>(game.vertexcount()); }

    template <typename F>
    void successors(int u, F f) const
    {
        if (!allowed[u].empty()) {
            for (int w : allowed[u]) f(w);
        } else {
            for (const int* e = game.outs(u); *e != -1; e++) f(*e);
        }
    }

    int levels() const { return static_cast<int>(odd.size()); }

    const pg::Game& game;
    const std::vector<std::vector<int>>& allowed;
    std::vector<int> level;     // priority index of each vertex
    std::vector<bool> odd;      // is the priority with this index odd?
};

// A graph given by successor lists and vertex priorities.
class Explicit final {
public:
    Explicit(const std::vector<std::vector<int>>& successor, const std::vector<int>& priority) : successor(successor)
    {
        toLevels(priority, level, odd);
    }

    int vertices() const { return static_cast<int>(successor.size()); }

    template <typename F>
    void successors(int u, F f) const
    {
        for (int w : successor[u]) f(w);
    }

    const std::vector<std::vector<int>>& successor;
    std::vector<int> level;
    std::vector<bool> odd;
};

// Does adding the edge v -> w close a cycle with an odd maximal priority?
// Search the pairs (u, m): u is reachable from w by a walk whose maximal
// priority has index m; the edge is bad if v is reached with an odd maximum.
// Searching walks instead of simple cycles is exact here: the restricted
// game has no cycle with an odd maximum yet, so if a closed walk through the
// new edge has an odd maximum, so does a simple cycle through it.
class CycleCheck final {
public:
    explicit CycleCheck(const Restricted& r)
        : r(r), words((r.levels() + 63) / 64), marks(static_cast<size_t>(r.game.vertexcount()) * words, 0)
    {
    }

    bool closesOddCycle(int v, int w)
    {
        bool found = false;
        visit(w, r.level[w]);
        while (!stack.empty() && !found) {
            const int u = stack.back().first;
            const int m = stack.back().second;
            stack.pop_back();
            if (u == v && r.odd[std::max(m, r.level[v])]) {
                found = true;
                break;
            }
            r.successors(u, [&](int x) { visit(x, std::max(m, r.level[x])); });
        }
        for (size_t i : touched) marks[i] = 0;
        touched.clear();
        stack.clear();
        return found;
    }

private:
    void visit(int u, int m)
    {
        const size_t i = static_cast<size_t>(u) * words + static_cast<size_t>(m / 64);
        const uint64_t bit = uint64_t(1) << (m % 64);
        if (marks[i] & bit) return;
        if (marks[i] == 0) touched.push_back(i);
        marks[i] |= bit;
        stack.emplace_back(u, m);
        work++;
    }

public:
    long work = 0;      // (vertex, level) pairs visited by all checks
private:

    const Restricted& r;
    const int words;
    std::vector<uint64_t> marks;
    std::vector<size_t> touched;
    std::vector<std::pair<int, int>> stack;
};

// The vertices reachable from the initial vertex.
template <typename Graph>
std::vector<int> reachable(const Graph& r, int initial)
{
    std::vector<bool> seen(r.vertices(), false);
    std::vector<int> result{initial};
    seen[initial] = true;
    for (size_t i = 0; i < result.size(); i++) {
        r.successors(result[i], [&](int x) {
            if (!seen[x]) {
                seen[x] = true;
                result.push_back(x);
            }
        });
    }
    return result;
}

// A cycle in the strongly connected component through one of its vertices
// with the given level (the maximal one): a shortest path back to it.
template <typename Graph>
std::vector<int> cycleThrough(const Graph& r, const std::vector<int>& component, int top)
{
    const int n = r.vertices();
    std::vector<char> inside(n, 0);
    for (int x : component) inside[x] = 1;
    int v = component[0];
    for (int x : component) {
        if (r.level[x] == top) v = x;
    }
    std::vector<int> parent(n, -1);
    std::deque<int> queue{v};
    while (!queue.empty()) {
        const int u = queue.front();
        queue.pop_front();
        int last = -1;
        r.successors(u, [&](int x) {
            if (!inside[x] || last != -1) return;
            if (x == v) {
                last = u;
            } else if (parent[x] == -1) {
                parent[x] = u;
                queue.push_back(x);
            }
        });
        if (last != -1) {
            std::vector<int> cycle;
            for (int x = last; x != v; x = parent[x]) cycle.push_back(x);
            cycle.push_back(v);
            std::reverse(cycle.begin(), cycle.end());
            return cycle;
        }
    }
    return {};
}

// Does the subgraph induced by the vertices in region have a cycle whose
// maximal priority is odd? Split it into strongly connected components; a
// component with an odd maximum has such a cycle, otherwise remove its
// vertices with the (even) maximum and examine what remains. If cycle is
// given, it receives such a cycle.
template <typename Graph>
bool hasOddCycle(const Graph& r, std::vector<int> region, std::vector<int>* cycle = nullptr)
{
    const int n = r.vertices();
    std::vector<int> inRegion(n, -1);     // region number, to restrict the graph
    std::vector<int> index(n, -1), low(n, 0);
    std::vector<bool> onStack(n, false);
    std::deque<std::vector<int>> work;
    work.push_back(std::move(region));
    int stamp = 0;
    while (!work.empty()) {
        std::vector<int> part = std::move(work.front());
        work.pop_front();
        stamp++;
        for (int v : part) {
            inRegion[v] = stamp;
            index[v] = -1;
        }
        // iterative Tarjan
        int counter = 0;
        std::vector<int> sccStack;
        std::vector<std::pair<int, std::vector<int>>> call;  // vertex, successors in the region
        for (int root : part) {
            if (index[root] != -1) continue;
            auto push = [&](int v) {
                index[v] = low[v] = counter++;
                sccStack.push_back(v);
                onStack[v] = true;
                std::vector<int> next;
                r.successors(v, [&](int x) {
                    if (inRegion[x] == stamp) next.push_back(x);
                });
                call.emplace_back(v, std::move(next));
            };
            push(root);
            while (!call.empty()) {
                const int v = call.back().first;
                std::vector<int>& next = call.back().second;
                if (!next.empty()) {
                    const int x = next.back();
                    next.pop_back();
                    if (index[x] == -1) {
                        push(x);
                    } else if (onStack[x]) {
                        low[v] = std::min(low[v], index[x]);
                    }
                    continue;
                }
                call.pop_back();
                if (!call.empty()) low[call.back().first] = std::min(low[call.back().first], low[v]);
                if (low[v] != index[v]) continue;
                // v is the root of a component
                std::vector<int> component;
                int x;
                do {
                    x = sccStack.back();
                    sccStack.pop_back();
                    onStack[x] = false;
                    component.push_back(x);
                } while (x != v);
                bool cyclic = component.size() > 1;
                if (!cyclic) {
                    r.successors(v, [&](int y) { cyclic = cyclic || y == v; });
                }
                if (!cyclic) continue;
                int top = 0;
                for (int y : component) top = std::max(top, r.level[y]);
                if (r.odd[top]) {
                    if (cycle != nullptr) *cycle = cycleThrough(r, component, top);
                    return true;
                }
                std::vector<int> rest;
                for (int y : component) {
                    if (r.level[y] != top) rest.push_back(y);
                }
                if (!rest.empty()) work.push_back(std::move(rest));
            }
        }
    }
    return false;
}

}  // namespace

void PermissiveStrategy::multiSolver(pg::Oink& engine, bool justifications)
{
    engine.setSolver(justifications ? "fpjm" : "fpim");
    engine.setRemoveLoops(false);
    engine.setRemoveWCWC(false);
    engine.setSolveSingle(false);
    engine.setBottomSCC(false);
}

std::vector<std::vector<int>> PermissiveStrategy::moves(const pg::Game& game, const pg::Solution& solution)
{
    std::vector<std::vector<int>> allowed(game.vertexcount());
    for (int v = 0; v < game.vertexcount(); v++) {
        if (game.owner(v) == 0 && solution.is_solved(v) && solution.winner(v) == 0) {
            pg::strategy_targets(game, solution, v, allowed[v]);
        }
    }
    return allowed;
}

std::vector<std::vector<int>> PermissiveStrategy::extend(const pg::Game& game, const pg::Solution& solution, int initial,
                                                         std::vector<std::vector<int>> allowed, long budget,
                                                         Stats& stats)
{
    const int n = game.vertexcount();
    Restricted r(game, allowed);
    CycleCheck check(r);

    std::vector<bool> seen(n, false);
    std::deque<int> queue{initial};
    seen[initial] = true;
    while (!queue.empty()) {
        const int v = queue.front();
        queue.pop_front();
        if (!allowed[v].empty()) {
            stats.vertices++;
            for (const int* e = game.outs(v); *e != -1; e++) {
                const int w = *e;
                if (solution.winner(w) != 0) continue;
                if (std::find(allowed[v].begin(), allowed[v].end(), w) != allowed[v].end()) continue;
                if (budget > 0 && check.work >= budget) {
                    stats.stopped = true;
                    break;
                }
                stats.candidates++;
                if (check.closesOddCycle(v, w)) continue;
                allowed[v].push_back(w);
                stats.added++;
            }
        }
        r.successors(v, [&](int x) {
            if (!seen[x]) {
                seen[x] = true;
                queue.push_back(x);
            }
        });
    }
    stats.work = check.work;
    return allowed;
}

bool PermissiveStrategy::verify(const pg::Game& game, const pg::Solution& solution, int initial,
                                const std::vector<std::vector<int>>& allowed)
{
    Restricted r(game, allowed);
    const std::vector<int> region = reachable(r, initial);
    for (int v : region) {
        if (solution.winner(v) != 0) return false;
    }
    return !hasOddCycle(r, region);
}

bool hasOddMaxCycle(const std::vector<std::vector<int>>& successors, const std::vector<int>& priorities, int initial)
{
    const Explicit graph(successors, priorities);
    return hasOddCycle(graph, reachable(graph, initial));
}

std::vector<int> oddMaxCycle(const std::vector<std::vector<int>>& successors, const std::vector<int>& priorities, int initial)
{
    const Explicit graph(successors, priorities);
    std::vector<int> cycle;
    hasOddCycle(graph, reachable(graph, initial), &cycle);
    return cycle;
}
