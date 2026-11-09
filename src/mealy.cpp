/**
 * Copyright Tom van Dijk
 */

#include <algorithm>
#include <cstdint>
#include <deque>
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <bddtools.hpp>
#include <mealy.hpp>
#include <permissive.hpp>
#include <statesearch.hpp>

#include <oink/game.hpp>
#include <oink/oink.hpp>

// ABC's headers have unused variables
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
extern "C" {
#include "sat/bsat/satSolver.h"
}
#pragma GCC diagnostic pop

using namespace sylvan;

namespace {

// BDDs that are protected from garbage collection until destruction.
class BddPool final {
public:
    BddPool() = default;
    BddPool(const BddPool&) = delete;
    BddPool& operator=(const BddPool&) = delete;

    ~BddPool()
    {
        for (MTBDD& b : items) mtbdd_unprotect(&b);
    }

    MTBDD keep(MTBDD b)
    {
        items.push_back(b);
        mtbdd_protect(&items.back());
        return b;
    }

private:
    std::deque<MTBDD> items;    // stable addresses
};

struct Option {
    int outputs;        // index of the set of allowed outputs
    int successor;      // index of the successor state
    int priority;       // priority of the move
};

// A nondeterministic Mealy machine: options[state][letter].
struct Machine {
    int initial = 0;
    std::vector<std::vector<std::vector<Option>>> options;

    int states() const { return static_cast<int>(options.size()); }
};

struct Aborted {
    const char* reason;
};

// Limits on the sizes for which the lower bound and the SAT encoding are tried,
// and on the work for the compatibility relation (pairs of options compared, a
// few seconds), so that the minimisation does not depend on time.
constexpr int MAX_COMPAT_STATES = 1500;
constexpr long MAX_COMPAT_WORK = 1000000000;
constexpr long MAX_CLAUSES = 20000000;
// Monitor codes are tried if (inputs x controller states x monitors) is at
// most MAX_MONITOR_WORK (a few seconds); splitting controller states for
// their flags may add up to a share of 1/MONITOR_SPLIT_SHARE of the states
// (each added state costs little, a flag can save much logic).
constexpr long MAX_MONITOR_WORK = 400000000;
constexpr int MONITOR_SPLIT_SHARE = 32;

// The node below the state variables for the given state (MSB first).
MTBDD cofactorState(MTBDD node, uint32_t firstVar, int bits, int state)
{
    for (int i = 0; i < bits && !mtbdd_isleaf(node); i++) {
        if (mtbdd_getvar(node) != firstVar + static_cast<uint32_t>(i)) continue;
        node = ((state >> (bits - 1 - i)) & 1) ? mtbdd_gethigh(node) : mtbdd_getlow(node);
    }
    return node;
}

int decodeBits(const uint8_t* bits, int count)
{
    int value = 0;
    for (int i = 0; i < count; i++) value = (value << 1) | (bits[i] == 1 ? 1 : 0);
    return value;
}

// Builds the machine from the strategy of the game.
class Extraction final {
public:
    Extraction(SymGame& game, MTBDD strategies, BddPool& pool) : g(game), strategies(strategies), pool(pool) {}

    Machine machine;
    std::vector<MTBDD> letterSets;  // inputs of each letter (BDD over the uap variables)
    std::vector<MTBDD> outputSets;  // allowed outputs of each option (BDD over the cap variables)

    void run()
    {
        findStates();
        findLetters();
        findOptions();
    }

    // The input valuations of every letter (empty if there are too many).
    std::vector<std::vector<std::vector<uint8_t>>> inputs() const
    {
        if (uapCount > MAX_INPUT_BITS) return {};
        std::vector<std::vector<std::vector<uint8_t>>> result(letterSets.size());
        std::vector<uint8_t> bits(uapCount);
        for (size_t l = 0; l < letterSets.size(); l++) {
            for (MTBDD lf = mtbdd_enum_all_first(letterSets[l], g.uap_vars, bits.data(), nullptr); lf != mtbdd_false;
                 lf = mtbdd_enum_all_next(letterSets[l], g.uap_vars, bits.data(), nullptr)) {
                result[l].push_back(bits);
            }
        }
        return result;
    }

private:
    void findStates()
    {
        // the moves allowed by the strategy, from state to next state
        MTBDD vars = pool.keep(mtbdd_set_addall(g.uap_vars, g.cap_vars));
        vars = pool.keep(mtbdd_set_addall(vars, g.np_vars));
        const MTBDD moves = pool.keep(sylvan_and_exists(strategies, g.trans, vars));
        const uint32_t firstS = mtbdd_set_first(g.s_vars);
        index[0] = 0;
        numbers.push_back(0);
        std::vector<uint8_t> bits(g.statebits);
        for (size_t i = 0; i < numbers.size(); i++) {
            const MTBDD next = cofactorState(moves, firstS, g.statebits, numbers[i]);
            for (MTBDD lf = mtbdd_enum_all_first(next, g.ns_vars, bits.data(), nullptr); lf != mtbdd_false;
                 lf = mtbdd_enum_all_next(next, g.ns_vars, bits.data(), nullptr)) {
                const int n = decodeBits(bits.data(), g.statebits);
                if (index.emplace(n, static_cast<int>(numbers.size())).second) numbers.push_back(n);
            }
        }
        for (int n : numbers) {
            tops.push_back(cofactorState(g.trans, firstS, g.statebits, n));
            tops.push_back(cofactorState(strategies, firstS, g.statebits, n));
        }
    }

    // The letters are the classes of inputs that every state treats alike:
    // follow the uap variables in all states at once, level by level, with
    // the inputs that reach each tuple of nodes (as a BDD), so that the work
    // is in the number of distinct tuples, not of paths. The letters are
    // numbered by their least input valuation (the order in which a
    // depth-first search, low first, meets them).
    void findLetters()
    {
        uapCount = static_cast<size_t>(g.uap_count);
        if (uapCount > 0) firstUap = mtbdd_set_first(g.uap_vars);
        BddPool local;                      // the inputs of the tuples of the levels
        std::map<std::vector<MTBDD>, MTBDD> level{{tops, mtbdd_true}};
        MTBDD literal = mtbdd_false, part = mtbdd_false;
        mtbdd_refs_pushptr(&literal);
        mtbdd_refs_pushptr(&part);
        auto merge = [&](std::map<std::vector<MTBDD>, MTBDD>& into, const std::vector<MTBDD>& nodes, MTBDD inputs) {
            auto it = into.find(nodes);
            if (it == into.end()) {
                into.emplace(nodes, local.keep(inputs));
            } else {
                it->second = local.keep(sylvan_or(it->second, inputs));
            }
        };
        for (size_t i = 0; i < uapCount; i++) {
            const uint32_t var = firstUap + static_cast<uint32_t>(i);
            std::map<std::vector<MTBDD>, MTBDD> next;
            for (const auto& [nodes, inputs] : level) {
                bool depends = false;
                for (MTBDD n : nodes) depends = depends || (!mtbdd_isleaf(n) && mtbdd_getvar(n) == var);
                if (!depends) {
                    merge(next, nodes, inputs);
                    continue;
                }
                std::vector<MTBDD> low(nodes), high(nodes);
                for (size_t k = 0; k < nodes.size(); k++) {
                    if (!mtbdd_isleaf(nodes[k]) && mtbdd_getvar(nodes[k]) == var) {
                        low[k] = mtbdd_getlow(nodes[k]);
                        high[k] = mtbdd_gethigh(nodes[k]);
                    }
                }
                literal = sylvan_nithvar(var);
                part = sylvan_and(inputs, literal);
                merge(next, low, part);
                literal = sylvan_ithvar(var);
                part = sylvan_and(inputs, literal);
                merge(next, high, part);
            }
            if (static_cast<long>(next.size()) > MAX_CLASSES) {
                mtbdd_refs_popptr(2);
                throw Aborted{"too many classes of inputs"};
            }
            level = std::move(next);
        }
        mtbdd_refs_popptr(2);
        // number the letters by their least input valuation
        std::vector<std::pair<std::vector<uint8_t>, const std::pair<const std::vector<MTBDD>, MTBDD>*>> letters;
        std::vector<uint8_t> bits(uapCount);
        for (const auto& entry : level) {
            mtbdd_enum_all_first(entry.second, g.uap_vars, bits.data(), nullptr);
            letters.emplace_back(bits, &entry);
        }
        std::sort(letters.begin(), letters.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
        for (const auto& [least, entry] : letters) {
            letterNodes.push_back(entry->first);
            letterSets.push_back(pool.keep(entry->second));
        }
    }

    // The options of a state on a letter: for every next priostate that
    // the allowed outputs reach, the outputs that reach it.
    void findOptions()
    {
        const uint32_t firstNp = mtbdd_set_first(g.np_vars);
        const int states = static_cast<int>(numbers.size());
        machine.options.assign(states, std::vector<std::vector<Option>>(letterNodes.size()));
        std::map<std::pair<MTBDD, MTBDD>, std::vector<Option>> cache;
        std::vector<uint8_t> bits(g.priobits + g.statebits);
        MTBDD reach = mtbdd_false;
        MTBDD outputs = mtbdd_false;
        mtbdd_refs_pushptr(&reach);
        mtbdd_refs_pushptr(&outputs);
        for (size_t l = 0; l < letterNodes.size(); l++) {
            for (int s = 0; s < states; s++) {
                const MTBDD after = letterNodes[l][2 * s];      // trans after the inputs
                const MTBDD allowed = letterNodes[l][2 * s + 1]; // allowed outputs
                auto it = cache.find({after, allowed});
                if (it == cache.end()) {
                    std::vector<Option> options;
                    for (MTBDD target : BDDTools::collectSubroots(after, firstNp)) {
                        if (target == mtbdd_false) continue;
                        reach = BDDTools::pathsToSubroot(after, firstNp, target);
                        outputs = sylvan_and(reach, allowed);
                        if (outputs == mtbdd_false) continue;
                        mtbdd_enum_all_first(target, g.pns_vars, bits.data(), nullptr);
                        const int priority = decodeBits(bits.data(), g.priobits);
                        const int next = decodeBits(bits.data() + g.priobits, g.statebits);
                        auto succ = index.find(next);
                        if (succ == index.end()) {
                            mtbdd_refs_popptr(2);
                            throw Aborted{"a successor of the strategy is not reachable"};
                        }
                        options.push_back(Option{outputId(outputs), succ->second, priority});
                    }
                    it = cache.emplace(std::make_pair(after, allowed), std::move(options)).first;
                }
                if (it->second.empty()) {
                    mtbdd_refs_popptr(2);
                    throw Aborted{"the strategy has no move for some input"};
                }
                machine.options[s][l] = it->second;
            }
        }
        mtbdd_refs_popptr(2);
    }

    int outputId(MTBDD outputs)
    {
        const auto inserted = outputIndex.emplace(outputs, static_cast<int>(outputSets.size()));
        if (inserted.second) outputSets.push_back(pool.keep(outputs));
        return inserted.first->second;
    }

    static constexpr long MAX_CLASSES = long(1) << 20;
    static constexpr size_t MAX_INPUT_BITS = 14;

    SymGame& g;
    const MTBDD strategies;             // the allowed moves (protected by the caller)
    BddPool& pool;
    std::map<int, int> index;           // state number -> state
    std::vector<int> numbers;           // state -> state number
    std::vector<MTBDD> tops;            // per state: trans and strategies below the state variables
    size_t uapCount = 0;
    uint32_t firstUap = 0;
    std::vector<std::vector<MTBDD>> letterNodes;
    std::map<MTBDD, int> outputIndex;
};

// Is the machine a winning multi-strategy, i.e., does every play that only
// makes allowed moves have an even maximal priority? (Moves are vertices with
// their priority between the states.)
bool winning(const Machine& m)
{
    std::vector<std::vector<int>> successors(m.states());
    std::vector<int> priorities(m.states(), 0);
    std::map<std::pair<std::pair<int, int>, int>, int> moves;
    for (int q = 0; q < m.states(); q++) {
        for (const std::vector<Option>& options : m.options[q]) {
            for (const Option& o : options) {
                const auto key = std::make_pair(std::make_pair(q, o.successor), o.priority);
                const auto inserted = moves.emplace(key, static_cast<int>(successors.size()));
                if (!inserted.second) continue;
                successors[q].push_back(inserted.first->second);
                successors.push_back({o.successor});
                priorities.push_back(o.priority);
            }
        }
    }
    return !hasOddMaxCycle(successors, priorities, m.initial);
}

// Merge bisimilar states: states with the same options (outputs and class of
// the successor, and the priority if it matters) on every letter. This is
// exact: a controller state can serve a state iff it can serve every state
// bisimilar to it.
Machine quotient(const Machine& m, std::vector<int>& block, bool priorities)
{
    const int n = m.states();
    block.assign(n, 0);
    int blocks = 1;
    for (;;) {
        std::map<std::vector<int>, int> ids;
        std::vector<int> next(n);
        for (int s = 0; s < n; s++) {
            std::vector<int> sig{block[s]};
            for (const std::vector<Option>& options : m.options[s]) {
                std::vector<std::tuple<int, int, int>> moves;
                for (const Option& o : options) moves.emplace_back(o.outputs, block[o.successor], priorities ? o.priority : 0);
                std::sort(moves.begin(), moves.end());
                sig.push_back(-1);
                for (const auto& [outputs, successor, priority] : moves) {
                    sig.push_back(outputs);
                    sig.push_back(successor);
                    if (priorities) sig.push_back(priority);
                }
            }
            next[s] = ids.emplace(std::move(sig), static_cast<int>(ids.size())).first->second;
        }
        block = std::move(next);
        if (static_cast<int>(ids.size()) == blocks) break;
        blocks = static_cast<int>(ids.size());
    }
    Machine q;
    q.initial = block[m.initial];
    q.options.resize(blocks);
    std::vector<bool> done(blocks, false);
    for (int s = 0; s < n; s++) {
        if (done[block[s]]) continue;
        done[block[s]] = true;
        q.options[block[s]] = m.options[s];
        for (std::vector<Option>& options : q.options[block[s]]) {
            for (Option& o : options) o.successor = block[o.successor];
        }
    }
    return q;
}

// Do two output sets overlap?
class Overlap final {
public:
    explicit Overlap(const std::vector<MTBDD>& sets) : sets(sets) {}

    bool operator()(int a, int b)
    {
        if (a == b) return true;
        const auto key = std::make_pair(std::min(a, b), std::max(a, b));
        auto it = cache.find(key);
        if (it != cache.end()) return it->second;
        const bool result = sylvan_and(sets[a], sets[b]) != mtbdd_false;
        cache.emplace(key, result);
        return result;
    }

private:
    const std::vector<MTBDD>& sets;
    std::map<std::pair<int, int>, bool> cache;
};

// The states that every controller visits: q is unavoidable if the
// environment can force a visit to q whatever allowed outputs the controller
// chooses. reach[s], the states that the environment can force from s, is the
// least set with s itself and, for some letter, the states that it can force
// from the successors of all options of s (computed for all states at once,
// as bitsets). Stops after MAX_FORCED_WORK words; every state found so far is
// unavoidable. Without strong, only letters on which all options have the same
// successor count (the states forced along single successors).
constexpr long MAX_FORCED_WORK = 400000000;

std::vector<bool> forcedStates(const Machine& m, bool strong)
{
    const int n = m.states();
    const int words = (n + 63) / 64;
    std::vector<uint64_t> reach(static_cast<size_t>(n) * words, 0);
    auto row = [&](int q) { return reach.data() + static_cast<size_t>(q) * words; };
    for (int q = 0; q < n; q++) row(q)[q / 64] |= uint64_t(1) << (q % 64);
    std::vector<uint64_t> meet(words);
    long work = 0;
    for (bool changed = true; changed && work < MAX_FORCED_WORK;) {
        changed = false;
        for (int q = 0; q < n; q++) {
            uint64_t* target = row(q);
            for (const std::vector<Option>& options : m.options[q]) {
                if (!strong && std::any_of(options.begin(), options.end(),
                                           [&](const Option& o) { return o.successor != options[0].successor; })) {
                    continue;
                }
                const uint64_t* first = row(options[0].successor);
                std::copy(first, first + words, meet.begin());
                for (size_t o = 1; o < options.size(); o++) {
                    const uint64_t* other = row(options[o].successor);
                    for (int w = 0; w < words; w++) meet[w] &= other[w];
                }
                for (int w = 0; w < words; w++) {
                    if (meet[w] & ~target[w]) {
                        target[w] |= meet[w];
                        changed = true;
                    }
                }
                work += static_cast<long>(options.size() + 1) * words;
            }
        }
    }
    std::vector<bool> forced(n);
    const uint64_t* initial = row(m.initial);
    for (int q = 0; q < n; q++) forced[q] = (initial[q / 64] >> (q % 64)) & 1;
    return forced;
}

// Adds the states that every winning controller visits, for moves that are
// not all winning (the search over all winning moves): q is unavoidable if
// the controller loses the parity game of the machine in which the moves
// into q lose. One game per state (solved by Oink), except for the states
// known to be unavoidable and those that the winning strategy of an earlier
// game avoids. Stops before the edges of the games would exceed
// MAX_WINNING_WORK (a bound on their size, not on Oink's work for solving them);
// the states not decided by then count as avoidable, which is sound.
constexpr long MAX_WINNING_WORK = 100000000;

void addWinningUnavoidable(const Machine& m, std::vector<bool>& unavoidable)
{
    // vertices: the states (the environment chooses a letter), per state and
    // letter one of the controller (it chooses an option), per priority and
    // successor a move, and a losing sink for the moves into q
    const int n = m.states();
    int count = n;
    std::map<std::pair<int, int>, int> moves;      // (priority, successor) -> vertex
    for (const std::vector<std::vector<Option>>& letters : m.options) {
        count += static_cast<int>(letters.size());
        for (const std::vector<Option>& options : letters) {
            for (const Option& o : options) moves.emplace(std::make_pair(o.priority, o.successor), 0);
        }
    }
    int maxPriority = 0;
    for (auto& [move, vertex] : moves) {
        vertex = count++;
        maxPriority = std::max(maxPriority, move.first);
    }
    const int sink = count++;
    std::vector<std::vector<std::vector<int>>> choices(n);    // per state and letter: the moves
    long edges = static_cast<long>(moves.size()) + 1;
    for (int s = 0; s < n; s++) {
        for (const std::vector<Option>& options : m.options[s]) {
            std::vector<int> targets;
            for (const Option& o : options) targets.push_back(moves.at(std::make_pair(o.priority, o.successor)));
            std::sort(targets.begin(), targets.end());
            targets.erase(std::unique(targets.begin(), targets.end()), targets.end());
            edges += static_cast<long>(targets.size()) + 1;
            choices[s].push_back(std::move(targets));
        }
    }

    std::vector<bool> avoidable(n, false);
    std::vector<int> mapping(count);
    std::ostream discard(nullptr);
    long work = 0;
    for (int q = 0; q < n && work + edges <= MAX_WINNING_WORK; q++) {
        if (unavoidable[q] || avoidable[q]) continue;
        work += edges;
        pg::Game game(count);
        int v = n;
        for (int s = 0; s < n; s++) {
            game.init_vertex(s, 0, 1);
            game.e_start(s);
            for (size_t l = 0; l < choices[s].size(); l++) game.e_add(s, v + static_cast<int>(l));
            game.e_finish();
            for (const std::vector<int>& targets : choices[s]) {
                game.init_vertex(v, 0, 0);
                game.e_start(v);
                for (int t : targets) game.e_add(v, t);
                game.e_finish();
                v++;
            }
        }
        for (const auto& [move, vertex] : moves) {
            game.init_vertex(vertex, move.first, 1);
            game.e_start(vertex);
            game.e_add(vertex, move.second == q ? sink : move.second);
            game.e_finish();
        }
        game.init_vertex(sink, maxPriority % 2 == 1 ? maxPriority : maxPriority + 1, 1);
        game.e_start(sink);
        game.e_add(sink, sink);
        game.e_finish();
        // Oink orders the vertices by priority; keep the numbering, of the game and its solution
        game.sort(mapping.data());
        pg::Oink engine(game, discard);
        engine.setSolver("tl");
        engine.run();
        pg::permute(game, engine.solution(), mapping.data());
        const pg::Solution& solution = engine.solution();
        if (solution.winner(m.initial) != 0) {
            unavoidable[q] = true;
            continue;
        }
        // the winning strategy avoids q and every state that it does not reach
        std::vector<char> seen(count, 0);
        std::vector<int> todo{m.initial};
        seen[m.initial] = 1;
        while (!todo.empty()) {
            const int u = todo.back();
            todo.pop_back();
            const int chosen = game.owner(u) == 0 ? solution.strategy(u) : -1;
            for (const int* e = game.outs(u); *e != -1; e++) {
                if (chosen != -1 && *e != chosen) continue;
                if (!seen[*e]) {
                    seen[*e] = 1;
                    todo.push_back(*e);
                }
            }
        }
        for (int s = 0; s < n; s++) {
            if (!seen[s]) avoidable[s] = true;
        }
    }
}

// A set of pairwise incompatible states that every controller visits, so
// they need distinct controller states: two states are compatible if, on
// every letter, they have options with overlapping outputs whose successors
// are compatible (greatest fixed point). With winners (for moves that are not
// all winning), the states that every winning controller visits. Also gives
// the compatibility relation, or an over-approximation of it if this takes
// more than MAX_COMPAT_WORK (then the set is just the initial state), or
// nothing for large machines.
std::vector<int> incompatibleStates(const Machine& m, Overlap& overlap, bool strong, bool winners,
                                    std::vector<std::vector<char>>& compatible)
{
    const int n = m.states();
    std::vector<int> clique{m.initial};
    if (n > MAX_COMPAT_STATES) return clique;
    std::vector<bool> forced = forcedStates(m, strong);
    if (winners) addWinningUnavoidable(m, forced);
    compatible.assign(n, std::vector<char>(n, 1));
    long work = 0;
    for (bool changed = true; changed;) {
        changed = false;
        for (int a = 0; a < n; a++) {
            if (work > MAX_COMPAT_WORK) return clique;
            for (int b = a + 1; b < n; b++) {
                if (!compatible[a][b]) continue;
                bool ok = true;
                for (size_t l = 0; l < m.options[a].size() && ok; l++) {
                    work += static_cast<long>(m.options[a][l].size() * m.options[b][l].size());
                    ok = false;
                    for (const Option& oa : m.options[a][l]) {
                        for (const Option& ob : m.options[b][l]) {
                            if (compatible[oa.successor][ob.successor] && overlap(oa.outputs, ob.outputs)) {
                                ok = true;
                                break;
                            }
                        }
                        if (ok) break;
                    }
                }
                if (!ok) {
                    compatible[a][b] = compatible[b][a] = 0;
                    changed = true;
                }
            }
        }
    }
    // greedily, the states incompatible with most of the others first; the clique need not contain
    // the initial state (the initial state is often compatible with all others), but it is the first
    // anchor if that gives as large a clique
    auto greedy = [&](std::vector<int> candidates, std::vector<int> result) {
        while (!candidates.empty()) {
            int best = candidates[0];
            int bestDegree = -1;
            for (int c : candidates) {
                int degree = 0;
                for (int d : candidates) degree += compatible[c][d] ? 0 : 1;
                if (degree > bestDegree) {
                    best = c;
                    bestDegree = degree;
                }
            }
            result.push_back(best);
            std::vector<int> rest;
            for (int c : candidates) {
                if (c != best && !compatible[c][best]) rest.push_back(c);
            }
            candidates = std::move(rest);
        }
        return result;
    };
    std::vector<int> withInitial, all;
    for (int s = 0; s < n; s++) {
        if (!forced[s]) continue;
        all.push_back(s);
        if (!compatible[s][m.initial]) withInitial.push_back(s);
    }
    std::vector<int> first = greedy(withInitial, clique);
    if (!strong) return first;
    std::vector<int> other = greedy(all, {});
    return other.size() > first.size() ? other : first;
}

// The deterministic controller: successor and output bits per state and letter.
struct Controller {
    int initial = 0;
    std::vector<std::vector<int>> successor;
    std::vector<std::vector<std::vector<uint8_t>>> output;
};

// The controller that the bisimulation quotient gives directly: in every
// state and on every letter, the first option, with its first output.
Controller fromQuotient(const Machine& q, const std::vector<MTBDD>& outputSets, SymGame& g)
{
    Controller ctl;
    ctl.initial = q.initial;
    const size_t letters = q.options[0].size();
    ctl.successor.assign(q.states(), std::vector<int>(letters, 0));
    ctl.output.assign(q.states(), std::vector<std::vector<uint8_t>>(letters, std::vector<uint8_t>(g.cap_count, 0)));
    for (int b = 0; b < q.states(); b++) {
        for (size_t l = 0; l < letters; l++) {
            const Option& o = q.options[b][l][0];
            ctl.successor[b][l] = o.successor;
            if (g.cap_count == 0) continue;
            std::vector<uint8_t>& bits = ctl.output[b][l];
            mtbdd_enum_first(outputSets[o.outputs], g.cap_vars, bits.data(), nullptr);
            for (uint8_t& bit : bits) bit = bit == 1 ? 1 : 0;
        }
    }
    return ctl;
}

// A move in the product of a controller with a machine: controller state c
// serves (quotient) state q and takes an option with the priority.
struct Step {
    int c, q, priority;
};

// Is the output (bits for the cap variables) in the set?
bool contains(MTBDD set, const std::vector<uint8_t>& bits, uint32_t firstCap)
{
    while (!mtbdd_isleaf(set)) {
        set = bits[mtbdd_getvar(set) - firstCap] ? mtbdd_gethigh(set) : mtbdd_getlow(set);
    }
    return set == mtbdd_true;
}

// A cycle with an odd maximal priority in the product of the controller with
// the machine, as its moves, or nothing if every play of the product is won.
std::vector<Step> losingCycle(const Machine& m, const std::vector<int>& block, const std::vector<MTBDD>& outputSets,
                              const Controller& ctl, uint32_t firstCap)
{
    std::map<std::pair<int, int>, int> pairs;
    std::vector<std::vector<int>> successors;
    std::vector<int> priorities;
    std::vector<Step> stepOf;                   // per vertex; c is -1 for pairs
    std::vector<std::pair<int, int>> todo;
    auto vertex = [&](int c, int q) {
        const auto inserted = pairs.emplace(std::make_pair(c, q), static_cast<int>(successors.size()));
        if (inserted.second) {
            successors.emplace_back();
            priorities.push_back(0);
            stepOf.push_back(Step{-1, 0, 0});
            todo.emplace_back(c, q);
        }
        return inserted.first->second;
    };
    vertex(ctl.initial, m.initial);
    while (!todo.empty()) {
        const auto [c, q] = todo.back();
        todo.pop_back();
        const int from = pairs.at({c, q});
        for (size_t l = 0; l < m.options[q].size(); l++) {
            const Option* move = nullptr;
            for (const Option& o : m.options[q][l]) {
                if (contains(outputSets[o.outputs], ctl.output[c][l], firstCap)) move = &o;
            }
            if (move == nullptr) throw Aborted{"the controller makes a move that the game does not allow"};
            const int mv = static_cast<int>(successors.size());
            successors.emplace_back();
            priorities.push_back(move->priority);
            stepOf.push_back(Step{c, block[q], move->priority});
            successors[from].push_back(mv);
            const int to = vertex(ctl.successor[c][l], move->successor);
            successors[mv].push_back(to);
        }
    }
    std::vector<Step> cycle;
    for (int v : oddMaxCycle(successors, priorities, 0)) {
        if (stepOf[v].c != -1) cycle.push_back(stepOf[v]);
    }
    return cycle;
}

// The SAT encoding for a controller with k states. The first ones are
// anchored: controller state i serves anchors[i], one of the pairwise
// incompatible states of the lower bound. A controller state that serves an
// anchor can only serve states compatible with it, which leaves most states
// few possible controller states (often one); the other controller states
// are free. Variables, only for possible combinations: x(q,c) "controller
// state c serves q", s(c,l,d) "c goes to d on letter l" (exactly one),
// y(c,l,j) "output bit j of c on l", n(c,l,q) "the successor of c on l must
// serve q", and Tseitin variables for "the output of c on l is in output set
// o". Letters on which all states that c can serve have the same options are
// interchangeable for c, so c has these variables per class of such letters.
class Encoding final {
public:
    Encoding(const Machine& m, const std::vector<MTBDD>& outputSets, int capCount, uint32_t firstCap,
             const std::vector<int>& anchors, const std::vector<std::vector<char>>& compatible, int k, bool priorities)
        : m(m), outputSets(outputSets), capCount(capCount), firstCap(firstCap), anchors(anchors), k(k),
          priorities(priorities),
          anchored(static_cast<int>(anchors.size())), states(m.states()),
          letters(static_cast<int>(m.options[0].size())), solver(sat_solver_new()),
          xVar(static_cast<size_t>(states) * k, -1), sVar(k), yStart(k, 0)
    {
        anchorsOf.resize(states);
        servable.resize(anchored);
        for (int q = 0; q < states; q++) {
            allStates.push_back(q);
            for (int i = 0; i < anchored; i++) {
                if (compatible.empty() || compatible[q][anchors[i]]) {
                    anchorsOf[q].push_back(i);
                    servable[i].push_back(q);
                }
            }
        }
        anySuccessor.resize(letters);
        for (int l = 0; l < letters; l++) anySuccessor[l] = successors(allStates, l);
        for (int i = 0; i < anchored; i++) anchorClasses.push_back(classify(servable[i]));
        if (k > anchored) freeClasses = classify(allStates);
    }

    Encoding(const Encoding&) = delete;
    Encoding& operator=(const Encoding&) = delete;
    ~Encoding() { sat_solver_delete(solver); }

    int variables() const { return next; }
    int clauses() const { return sat_solver_nclauses(solver); }

    // The number of clauses build() adds, apart from the Tseitin encoding of
    // the output sets.
    long estimatedSize() const
    {
        long size = 0;
        for (int c = 0; c < k; c++) {
            for (int l : classes(c).letter) {
                for (int q : served(c)) size += 1 + static_cast<long>(m.options[q][l].size());
                const std::vector<int>& need = c < anchored ? successors(servable[c], l) : anySuccessor[l];
                size += 1 + 3 * static_cast<long>(targets(need).size());
                for (int q : need) size += 1 + 2 * static_cast<long>(anchorsOf[q].size() + (k - anchored));
            }
        }
        return size;
    }

    // Add the clauses: l_False if they are trivially unsatisfiable, l_Undef
    // (and stops) once more than limit clauses are added, else l_True.
    int build(long limit)
    {
        trueVar = next++;
        for (int c = 0; c < k; c++) {
            yStart[c] = next;
            next += static_cast<int>(classes(c).letter.size()) * capCount;
        }
        bool ok = clause({pos(trueVar)});
        for (int i = 0; i < anchored && ok; i++) ok = clause({pos(x(anchors[i], i))});
        // the controller state that starts, serving the initial state: its anchor, or else one of the
        // anchors compatible with it or the first free state (the free states are interchangeable)
        initialAnchor = static_cast<int>(std::find(anchors.begin(), anchors.end(), m.initial) - anchors.begin());
        if (initialAnchor == anchored) {
            initialAnchor = -1;
            std::vector<lit> some;
            for (int d : anchorsOf[m.initial]) initialVar.emplace_back(d, next++);
            if (k > anchored) initialVar.emplace_back(anchored, next++);
            for (const auto& [d, v] : initialVar) {
                some.push_back(pos(v));
                ok = ok && clause({neg(v), pos(x(m.initial, d))});
            }
            ok = ok && !some.empty() && clause(some) && atMostOne(some);
        }
        for (int c = 0; c < k && ok; c++) {
            if (added > limit) return l_Undef;
            const std::vector<int>& letterOf = classes(c).letter;
            sVar[c].resize(letterOf.size());
            for (int cls = 0; cls < static_cast<int>(letterOf.size()) && ok; cls++) {
                const int l = letterOf[cls];
                const std::vector<int> need = c < anchored ? successors(servable[c], l) : anySuccessor[l];
                // exactly one successor, among the controller states that can serve a needed state
                std::vector<std::pair<int, int>>& succ = sVar[c][cls];
                for (int d : targets(need)) succ.emplace_back(d, next++);
                std::vector<lit> some;
                for (const auto& [d, v] : succ) some.push_back(pos(v));
                ok = !some.empty() && clause(some) && atMostOne(some);
                // if the successor must serve q, it is a controller state that does
                for (int q : need) {
                    const int nq = n(c, cls, q);
                    std::vector<lit> lits{neg(nq)};
                    for (const auto& [d, v] : succ) {
                        if (!possible(q, d)) continue;
                        lits.push_back(pos(v));
                        ok = ok && clause({neg(nq), neg(v), pos(x(q, d))});
                    }
                    ok = ok && clause(lits);
                }
                // if c serves q, its output is allowed and the successor serves the successor of q
                std::map<MTBDD, lit>& memo = memos[c * letters + cls];
                for (int q : served(c)) {
                    if (!ok) break;
                    const int xq = x(q, c);
                    std::vector<lit> allowed{neg(xq)};
                    for (const Option& o : m.options[q][l]) {
                        const lit in = member(c, cls, outputSets[o.outputs], memo);
                        allowed.push_back(in);
                        ok = ok && clause({neg(xq), lit_neg(in), pos(n(c, cls, o.successor))});
                    }
                    ok = ok && clause(allowed);
                }
            }
        }
        return !ok ? l_False : added > limit ? l_Undef : l_True;
    }

    // l_True, l_False, or l_Undef if the work (with the clauses) reaches
    // limit first. The solver stops on propagations only, so it runs in chunks
    // of at most CHUNK propagations, and the work is checked in between.
    int solve(long limit)
    {
        sat_solver_setnvars(solver, next);
        while (work() < limit) {
            const long chunk = std::min(CHUNK, limit - work());
            const ABC_INT64_T props = solver->stats.propagations, inspects = solver->stats.inspects;
            const int r = sat_solver_solve(solver, nullptr, nullptr, 0, chunk, 0, 0);
            propagations += solver->stats.propagations - props;
            inspections += solver->stats.inspects - inspects;
            if (r != l_Undef) return r;
        }
        return l_Undef;
    }

    // The work so far, which follows the time closely: the clauses added, and
    // the propagations of solve() with a fifth of the clauses it inspects (a
    // propagation inspects 2 to 13 clauses). Every call of solve() propagates,
    // which bounds the rounds of exclude() and solve().
    long work() const { return added + propagations + inspections / 5; }

    // Exclude a losing cycle of the product, with its maximal (odd) priority
    // p: no cycle has maximal priority p iff there is a rank per pair (state,
    // controller state) that does not increase along moves with a priority up
    // to p and decreases along moves with priority p. Require this for all
    // moves that the controller can make from the pairs on the cycle, so that
    // every other cycle through these pairs is excluded too.
    bool exclude(const std::vector<Step>& cycle)
    {
        int p = 0;
        for (const Step& step : cycle) p = std::max(p, step.priority);
        bool ok = true;
        for (const Step& step : cycle) {
            if (ok && ranked.insert(std::make_tuple(step.q, step.c, p)).second) ok = rankPair(step.q, step.c, p);
        }
        return ok;
    }

    Controller controller() const
    {
        Controller ctl;
        ctl.initial = initialAnchor;
        for (const auto& [d, v] : initialVar) {
            if (sat_solver_var_value(solver, v)) ctl.initial = d;
        }
        ctl.successor.assign(k, std::vector<int>(letters, 0));
        ctl.output.assign(k, std::vector<std::vector<uint8_t>>(letters, std::vector<uint8_t>(capCount, 0)));
        for (int c = 0; c < k; c++) {
            for (int l = 0; l < letters; l++) {
                const int cls = classes(c).of[l];
                for (const auto& [d, v] : sVar[c][cls]) {
                    if (sat_solver_var_value(solver, v)) ctl.successor[c][l] = d;
                }
                for (int j = 0; j < capCount; j++) ctl.output[c][l][j] = sat_solver_var_value(solver, y(c, cls, j)) ? 1 : 0;
            }
        }
        return ctl;
    }

private:
    static lit pos(int v) { return toLit(v); }
    static lit neg(int v) { return lit_neg(toLit(v)); }

    bool clause(std::vector<lit> lits)
    {
        added++;
        return sat_solver_addclause(solver, lits.data(), lits.data() + lits.size()) != 0;
    }

    // The classes of letters on which the given states have the same options.
    struct Classes {
        std::vector<int> of;        // per letter: its class
        std::vector<int> letter;    // per class: a letter in it
    };

    Classes classify(const std::vector<int>& from) const
    {
        Classes result;
        result.of.resize(letters);
        std::map<std::vector<int>, int> ids;
        std::vector<int> key;
        for (int l = 0; l < letters; l++) {
            key.clear();
            for (int q : from) {
                for (const Option& o : m.options[q][l]) {
                    key.push_back(o.outputs);
                    key.push_back(o.successor);
                    if (priorities) key.push_back(o.priority);
                }
                key.push_back(-1);
            }
            const auto inserted = ids.emplace(key, static_cast<int>(result.letter.size()));
            if (inserted.second) result.letter.push_back(l);
            result.of[l] = inserted.first->second;
        }
        return result;
    }

    // The states that controller state c can serve, and its classes of letters.
    const std::vector<int>& served(int c) const { return c < anchored ? servable[c] : allStates; }
    const Classes& classes(int c) const { return c < anchored ? anchorClasses[c] : freeClasses; }

    // Can controller state d serve q?
    bool possible(int q, int d) const
    {
        return d >= anchored || std::binary_search(anchorsOf[q].begin(), anchorsOf[q].end(), d);
    }

    // The successors of the given states on letter l.
    std::vector<int> successors(const std::vector<int>& from, int l) const
    {
        std::vector<int> result;
        for (int q : from) {
            for (const Option& o : m.options[q][l]) result.push_back(o.successor);
        }
        std::sort(result.begin(), result.end());
        result.erase(std::unique(result.begin(), result.end()), result.end());
        return result;
    }

    // The controller states that can serve one of the given states.
    std::vector<int> targets(const std::vector<int>& need) const
    {
        std::vector<int> result;
        for (int q : need) result.insert(result.end(), anchorsOf[q].begin(), anchorsOf[q].end());
        std::sort(result.begin(), result.end());
        result.erase(std::unique(result.begin(), result.end()), result.end());
        for (int d = anchored; d < k; d++) result.push_back(d);
        return result;
    }

    // The rank constraints for priority p on all moves from (q,c).
    bool rankPair(int q, int c, int p)
    {
        bool ok = true;
        const Classes& cl = classes(c);
        for (int cls = 0; cls < static_cast<int>(cl.letter.size()) && ok; cls++) {
            for (const Option& o : m.options[q][cl.letter[cls]]) {
                if (o.priority > p) continue;
                const lit in = member(c, cls, outputSets[o.outputs], memos[c * letters + cls]);
                for (const auto& [d, v] : sVar[c][cls]) {
                    if (ok && possible(o.successor, d)) ok = rankMove(q, c, in, v, o.successor, d, p, o.priority == p);
                }
            }
        }
        return ok;
    }

    // If c serves q and makes the move (output in, successor variable v), then
    // rank_p(q,c) >= rank_p(next,d), strictly if asked (unsigned, MSB first).
    bool rankMove(int q, int c, lit in, int v, int nextState, int d, int p, bool strict)
    {
        const lit made = pos(next++);
        bool ok = clause({neg(x(q, c)), neg(v), lit_neg(in), made});
        const std::vector<int> a = rank(q, c, p);
        const std::vector<int> b = rank(nextState, d, p);
        lit equal = made;               // made, and the bits so far are equal
        for (size_t i = 0; i < a.size() && ok; i++) {
            const lit further = pos(next++);
            ok = clause({lit_neg(equal), pos(a[i]), neg(b[i])})
                 && clause({lit_neg(equal), neg(a[i]), neg(b[i]), further})
                 && clause({lit_neg(equal), pos(a[i]), pos(b[i]), further});
            equal = further;
        }
        if (ok && strict) ok = clause({lit_neg(equal)});
        return ok;
    }

    // The rank variables of (q,c) for priority p (MSB first).
    std::vector<int> rank(int q, int c, int p)
    {
        std::vector<int>& bits = ranks[std::make_tuple(q, c, p)];
        if (bits.empty()) {
            int width = 1;
            while ((1L << width) <= static_cast<long>(states) * k) width++;
            for (int i = 0; i < width; i++) bits.push_back(next++);
        }
        return bits;
    }

    // Pairwise for few literals, otherwise a sequential counter.
    bool atMostOne(const std::vector<lit>& lits)
    {
        bool ok = true;
        if (lits.size() <= 5) {
            for (size_t i = 0; i < lits.size(); i++) {
                for (size_t j = i + 1; j < lits.size() && ok; j++) ok = clause({lit_neg(lits[i]), lit_neg(lits[j])});
            }
            return ok;
        }
        lit before = lit_neg(pos(trueVar));    // "one of the earlier literals holds"
        for (size_t i = 0; i < lits.size() && ok; i++) {
            const lit now = pos(next++);
            ok = clause({lit_neg(lits[i]), now}) && clause({lit_neg(before), now})
                 && clause({lit_neg(lits[i]), lit_neg(before)});
            before = now;
        }
        return ok;
    }

    int x(int q, int c)
    {
        int& v = xVar[static_cast<size_t>(q) * k + c];
        if (v == -1) v = next++;
        return v;
    }

    int n(int c, int cls, int q)
    {
        const auto inserted = nVar.emplace((static_cast<long>(c) * letters + cls) * states + q, next);
        if (inserted.second) next++;
        return inserted.first->second;
    }

    int y(int c, int cls, int j) const { return yStart[c] + cls * capCount + j; }

    // A literal that holds iff the output of c on the letters of class cls is in the set.
    lit member(int c, int cls, MTBDD set, std::map<MTBDD, lit>& memo)
    {
        if (set == mtbdd_true) return pos(trueVar);
        if (set == mtbdd_false) return neg(trueVar);
        auto it = memo.find(set);
        if (it != memo.end()) return it->second;
        const lit high = member(c, cls, mtbdd_gethigh(set), memo);
        const lit low = member(c, cls, mtbdd_getlow(set), memo);
        const lit bit = pos(y(c, cls, static_cast<int>(mtbdd_getvar(set) - firstCap)));
        const lit t = pos(next++);
        clause({lit_neg(bit), lit_neg(high), t});
        clause({lit_neg(bit), high, lit_neg(t)});
        clause({bit, lit_neg(low), t});
        clause({bit, low, lit_neg(t)});
        memo.emplace(set, t);
        return t;
    }

    const Machine& m;
    const std::vector<MTBDD>& outputSets;
    const int capCount;
    const uint32_t firstCap;
    const std::vector<int>& anchors;
    const int k;
    const bool priorities;                      // must letter classes agree on priorities?
    const int anchored;
    const int states;
    const int letters;
    sat_solver* solver;
    int initialAnchor = 0;                      // the controller state serving the initial state, or -1:
    std::vector<std::pair<int, int>> initialVar;// the candidates and their selector variables
    long added = 0;                             // clauses
    long propagations = 0;                      // of solve()
    long inspections = 0;                       // of clauses, by solve()
    static constexpr long CHUNK = 10000000;
    std::vector<std::vector<int>> anchorsOf;    // per state: the anchors it is compatible with
    std::vector<std::vector<int>> servable;     // per anchor: the states compatible with it
    std::vector<int> allStates;
    std::vector<std::vector<int>> anySuccessor; // per letter: the successors of all states
    std::vector<Classes> anchorClasses;         // per anchor
    Classes freeClasses;                        // for the free states
    int next = 0;
    int trueVar = 0;
    std::vector<int> xVar;                      // per (q,c), or -1
    std::vector<std::vector<std::vector<std::pair<int, int>>>> sVar;   // per c and class: the possible successors and their variables
    std::vector<int> yStart;                    // per c: the first output variable
    std::map<int, std::map<MTBDD, lit>> memos;  // per (c, class): the literals for output sets
    std::set<std::tuple<int, int, int>> ranked; // pairs (q,c) with rank constraints, per priority
    std::map<std::tuple<int, int, int>, std::vector<int>> ranks;
    std::unordered_map<long, int> nVar;
};

// The controller states reachable from the initial one, in breadth-first order.
std::vector<int> reachableStates(const Controller& ctl)
{
    std::vector<bool> seen(ctl.successor.size(), false);
    std::vector<int> order{ctl.initial};
    seen[ctl.initial] = true;
    for (size_t i = 0; i < order.size(); i++) {
        for (int d : ctl.successor[order[i]]) {
            if (!seen[d]) {
                seen[d] = true;
                order.push_back(d);
            }
        }
    }
    return order;
}

// Sets of outputs (BDDs over the cap variables), each stored once.
class OutputTable final {
public:
    explicit OutputTable(BddPool& pool) : pool(pool) {}

    int add(MTBDD set)
    {
        const auto inserted = index.emplace(set, static_cast<int>(items.size()));
        if (inserted.second) items.push_back(pool.keep(set));
        return inserted.first->second;
    }

    const std::vector<MTBDD>& sets() const { return items; }

private:
    BddPool& pool;
    std::vector<MTBDD> items;
    std::map<MTBDD, int> index;
};

// A move of the controller: any output in the set, then the successor.
struct Choice {
    int outputs;        // in the output table
    int successor;
};

// Per controller state and letter: the moves that the controller may make.
using Moves = std::vector<std::vector<std::vector<Choice>>>;

// The moves of the deterministic controller.
Moves movesOf(const Controller& ctl, OutputTable& table, MTBDD capVars)
{
    Moves moves(ctl.successor.size());
    for (size_t c = 0; c < ctl.successor.size(); c++) {
        for (size_t l = 0; l < ctl.successor[c].size(); l++) {
            std::vector<uint8_t> bits(ctl.output[c][l]);
            moves[c].push_back({Choice{table.add(sylvan_cube(capVars, bits.data())), ctl.successor[c][l]}});
        }
    }
    return moves;
}

// The pairs (controller state, state of the machine) that the controller can
// reach, whichever of its moves it makes; empty if it has no move for some
// letter in a pair it reaches, or one of its moves produces an output that the
// strategy does not allow.
std::vector<std::vector<char>> reached(const Machine& m, const std::vector<MTBDD>& optionSets,
                                       const std::vector<MTBDD>& moveSets, int initial, const Moves& moves)
{
    std::vector<std::vector<char>> seen(moves.size(), std::vector<char>(m.states(), 0));
    std::vector<std::pair<int, int>> todo{{initial, m.initial}};
    seen[initial][m.initial] = 1;
    MTBDD covered = mtbdd_false;
    mtbdd_refs_pushptr(&covered);
    bool ok = true;
    while (!todo.empty() && ok) {
        const auto [c, q] = todo.back();
        todo.pop_back();
        for (size_t l = 0; l < m.options[q].size() && ok; l++) {
            if (l >= moves[c].size() || moves[c][l].empty()) {
                ok = false;
                break;
            }
            for (const Choice& choice : moves[c][l]) {
                const MTBDD out = moveSets[choice.outputs];
                covered = mtbdd_false;
                for (const Option& o : m.options[q][l]) {
                    if (sylvan_and(out, optionSets[o.outputs]) == mtbdd_false) continue;
                    covered = sylvan_or(covered, optionSets[o.outputs]);
                    char& mark = seen[choice.successor][o.successor];
                    if (!mark) {
                        mark = 1;
                        todo.emplace_back(choice.successor, o.successor);
                    }
                }
                ok = ok && sylvan_and(out, sylvan_not(covered)) == mtbdd_false;
            }
        }
    }
    mtbdd_refs_popptr(1);
    if (!ok) seen.clear();
    return seen;
}

// Is every play of the product of the moves with the machine won (no cycle
// with an odd maximal priority), whichever of its moves the controller makes?
bool winningMoves(const Machine& m, const std::vector<MTBDD>& optionSets, const std::vector<MTBDD>& moveSets,
                  int initial, const Moves& moves)
{
    std::map<std::pair<int, int>, int> pairs;
    std::vector<std::vector<int>> successors;
    std::vector<int> priorities;
    std::vector<std::pair<int, int>> todo;
    auto vertex = [&](int c, int q) {
        const auto inserted = pairs.emplace(std::make_pair(c, q), static_cast<int>(successors.size()));
        if (inserted.second) {
            successors.emplace_back();
            priorities.push_back(0);
            todo.emplace_back(c, q);
        }
        return inserted.first->second;
    };
    vertex(initial, m.initial);
    while (!todo.empty()) {
        const auto [c, q] = todo.back();
        todo.pop_back();
        const int from = pairs.at({c, q});
        for (size_t l = 0; l < m.options[q].size(); l++) {
            for (const Choice& choice : moves[c][l]) {
                for (const Option& o : m.options[q][l]) {
                    if (sylvan_and(moveSets[choice.outputs], optionSets[o.outputs]) == mtbdd_false) continue;
                    const int mv = static_cast<int>(successors.size());
                    successors.emplace_back();
                    priorities.push_back(o.priority);
                    successors[from].push_back(mv);
                    const int to = vertex(choice.successor, o.successor);
                    successors[mv].push_back(to);
                }
            }
        }
    }
    return !hasOddMaxCycle(successors, priorities, 0);
}

// All moves that keep the controller within the given pairs: on letter l,
// controller state c may produce an output and go to d if, in every state
// that c serves, the output is in an option whose successor d serves. So
// whichever of these moves the controller makes, it only reaches these pairs,
// and the encoding may choose (see AIGEncoder::setDontCares).
Moves relaxedMoves(const Machine& m, const std::vector<MTBDD>& optionSets, size_t states,
                   const std::vector<std::vector<char>>& pairs, OutputTable& table)
{
    std::vector<std::vector<int>> served(states), servers(m.states());
    for (size_t c = 0; c < states; c++) {
        for (int q = 0; q < m.states(); q++) {
            if (!pairs[c][q]) continue;
            served[c].push_back(q);
            servers[q].push_back(static_cast<int>(c));
        }
    }
    const size_t letters = m.options[0].size();
    Moves moves(states);
    MTBDD all = mtbdd_false, some = mtbdd_false;
    mtbdd_refs_pushptr(&all);
    mtbdd_refs_pushptr(&some);
    for (size_t c = 0; c < states; c++) {
        if (served[c].empty()) continue;
        moves[c].resize(letters);
        for (size_t l = 0; l < letters; l++) {
            // the successor serves a successor of the first state that c serves
            std::vector<int> candidates;
            for (const Option& o : m.options[served[c][0]][l]) {
                candidates.insert(candidates.end(), servers[o.successor].begin(), servers[o.successor].end());
            }
            std::sort(candidates.begin(), candidates.end());
            candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
            for (int d : candidates) {
                all = mtbdd_true;
                for (size_t i = 0; i < served[c].size() && all != mtbdd_false; i++) {
                    some = mtbdd_false;
                    for (const Option& o : m.options[served[c][i]][l]) {
                        if (pairs[d][o.successor]) some = sylvan_or(some, optionSets[o.outputs]);
                    }
                    all = sylvan_and(all, some);
                }
                if (all != mtbdd_false) moves[c][l].push_back(Choice{table.add(all), d});
            }
        }
    }
    mtbdd_refs_popptr(2);
    return moves;
}

// The controller states reachable from the initial one, in breadth-first or
// depth-first order (following the letters in order).
std::vector<int> numbering(int initial, const Moves& moves, bool depthFirst)
{
    std::vector<char> seen(moves.size(), 0);
    std::vector<int> order;
    if (depthFirst) {
        std::vector<int> stack{initial};
        while (!stack.empty()) {
            const int c = stack.back();
            stack.pop_back();
            if (seen[c]) continue;
            seen[c] = 1;
            order.push_back(c);
            for (auto letter = moves[c].rbegin(); letter != moves[c].rend(); ++letter) {
                for (auto choice = letter->rbegin(); choice != letter->rend(); ++choice) {
                    if (!seen[choice->successor]) stack.push_back(choice->successor);
                }
            }
        }
        return order;
    }
    order.push_back(initial);
    seen[initial] = 1;
    for (size_t i = 0; i < order.size(); i++) {
        for (const std::vector<Choice>& choices : moves[order[i]]) {
            for (const Choice& choice : choices) {
                if (!seen[choice.successor]) {
                    seen[choice.successor] = 1;
                    order.push_back(choice.successor);
                }
            }
        }
    }
    return order;
}

// The codes of the controller states in the given order: their positions.
std::vector<int> positions(const std::vector<int>& order, size_t states, int& statebits)
{
    std::vector<int> codes(states, -1);
    for (size_t i = 0; i < order.size(); i++) codes[order[i]] = static_cast<int>(i);
    statebits = 1;
    while ((1 << statebits) < static_cast<int>(order.size())) statebits++;
    return codes;
}

// Codes as cubes: with n states and b bits, the 2^b - n states whose code
// with the most significant bit set is unused also get that code, as a cube
// with that bit free, so that no code is unused and the moves into these
// states leave that bit free. The states in breadth-first order get the
// codes in order, or with byIncoming, the states with the most moves into
// them get the cubes (the initial state keeps code 0). Returns the order of
// the states, with their codes, the number of bits and the free bits (by
// state), or nothing if every code is used.
std::vector<int> cubeNumbering(int initial, const Moves& moves, bool byIncoming, std::vector<int>& codes,
                               int& statebits, std::vector<uint32_t>& free)
{
    std::vector<int> order = numbering(initial, moves, false);
    codes = positions(order, moves.size(), statebits);
    const int n = static_cast<int>(order.size());
    const int half = 1 << (statebits - 1);
    if (n == 2 * half) return {};
    if (byIncoming) {
        std::vector<long> incoming(moves.size(), 0);
        for (const auto& letters : moves) {
            for (const auto& choices : letters) {
                for (const Choice& choice : choices) incoming[choice.successor]++;
            }
        }
        // the codes from n - half to half - 1 get cubes: give them to the states with the most
        // incoming moves (stable, so ties keep the breadth-first order), the others keep theirs in order
        std::vector<int> rest(order.begin() + 1, order.end());
        std::stable_sort(rest.begin(), rest.end(), [&](int a, int b) { return incoming[a] > incoming[b]; });
        const int cubes = 2 * half - n;
        std::vector<int> chosen(rest.begin(), rest.begin() + std::min<int>(cubes, static_cast<int>(rest.size())));
        std::vector<char> isChosen(moves.size(), 0);
        for (int c : chosen) isChosen[c] = 1;
        std::vector<int> others{initial};
        for (int c : order) {
            if (c != initial && !isChosen[c]) others.push_back(c);
        }
        // others take the codes below n - half and from half on, chosen the codes in between
        order.clear();
        size_t o = 0;
        for (int code = 0; code < n; code++) {
            const bool cube = code >= n - half && code < half;
            if (cube && !chosen.empty()) {
                order.push_back(chosen.front());
                chosen.erase(chosen.begin());
            } else {
                order.push_back(others[o++]);
            }
        }
        codes = positions(order, moves.size(), statebits);
    }
    free.assign(moves.size(), 0);
    for (int c : order) {
        if (codes[c] >= n - half && codes[c] < half) free[c] = static_cast<uint32_t>(half);
    }
    return order;
}

// Better state codes by adjacency (as in the classical heuristics for the state
// assignment of finite state machines): states that are related get codes
// that differ in few bits. Two states are related for every letter on which
// they can go to the same successor or give the same set of outputs (then their
// functions agree, and can share a cube of codes), and the successors of a
// state on different letters are related (then its next-state functions depend
// little on the input). Swaps the codes (statebits bits, by state) of two of
// the states in order (or of a state and an unused code) while that lowers the
// sum of the weights times the Hamming distances, at most MAX_ADJACENCY_PASSES
// times over all states; only the lowest free bits of a code may change (the
// others are flags), and the first state in order (the initial one) keeps its
// code; it stops after MAX_ADJACENCY_WORK terms of the cost (about half a
// second). False (and the codes unchanged) for more than MAX_ADJACENCY_STATES
// states.
constexpr int MAX_ADJACENCY_STATES = 1024;
constexpr int MAX_ADJACENCY_PASSES = 8;
constexpr long MAX_ADJACENCY_WORK = 200000000;

bool improveCodes(const Moves& moves, const std::vector<int>& order, std::vector<int>& codes, int statebits, int free)
{
    const int n = static_cast<int>(order.size());
    if (n > MAX_ADJACENCY_STATES || n < 3 || free < 1) return false;
    std::vector<int> index(moves.size(), -1);   // the position of each state in order
    for (int i = 0; i < n; i++) index[order[i]] = i;
    std::vector<std::vector<long>> weight(n, std::vector<long>(n, 0));
    auto relate = [&](int a, int b) {
        if (a == b || index[a] < 0 || index[b] < 0) return;
        weight[index[a]][index[b]]++;
        weight[index[b]][index[a]]++;
    };
    const size_t letters = moves[order[0]].size();
    for (size_t l = 0; l < letters; l++) {
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                bool successor = false, output = false;
                for (const Choice& x : moves[order[i]][l]) {
                    for (const Choice& y : moves[order[j]][l]) {
                        successor = successor || x.successor == y.successor;
                        output = output || x.outputs == y.outputs;
                    }
                }
                if (successor) relate(order[i], order[j]);
                if (output) relate(order[i], order[j]);
            }
        }
    }
    for (int i = 0; i < n; i++) {
        std::vector<int> next;
        for (size_t l = 0; l < letters; l++) {
            for (const Choice& x : moves[order[i]][l]) next.push_back(x.successor);
        }
        std::sort(next.begin(), next.end());
        next.erase(std::unique(next.begin(), next.end()), next.end());
        for (size_t a = 0; a < next.size(); a++) {
            for (size_t b = a + 1; b < next.size(); b++) relate(next[a], next[b]);
        }
    }
    // code[i]: the code of the state at position i; owner(c): the position with code c, or -1 (the
    // codes may have up to 30 bits, so the owners are a map of the used codes)
    std::vector<int> code(n);
    std::unordered_map<int, int> owners;
    for (int i = 0; i < n; i++) {
        code[i] = codes[order[i]];
        owners[code[i]] = i;
    }
    auto owner = [&](int c) {
        const auto it = owners.find(c);
        return it == owners.end() ? -1 : it->second;
    };
    auto distance = [](int x, int y) { return __builtin_popcount(static_cast<unsigned>(x ^ y)); };
    // the change of the cost if position i takes code c (and the owner of c, if any, takes i's code)
    auto delta = [&](int i, int c) {
        const int j = owner(c);
        long d = 0;
        for (int k = 0; k < n; k++) {
            if (k == i || k == j) continue;
            d += weight[i][k] * (distance(c, code[k]) - distance(code[i], code[k]));
            if (j >= 0) d += weight[j][k] * (distance(code[i], code[k]) - distance(c, code[k]));
        }
        return d;
    };
    const int span = 1 << free;
    long work = 0;
    for (int pass = 0; pass < MAX_ADJACENCY_PASSES && work < MAX_ADJACENCY_WORK; pass++) {
        bool improved = false;
        for (int i = 1; i < n && work < MAX_ADJACENCY_WORK; i++) {
            const int base = code[i] & ~(span - 1);
            int bestCode = -1;
            long bestDelta = 0;
            for (int c = base; c < base + span; c++) {
                if (c == code[i] || owner(c) == 0) continue;     // the first state keeps its code
                work += n;
                const long d = delta(i, c);
                if (d < bestDelta) {
                    bestDelta = d;
                    bestCode = c;
                }
            }
            if (bestCode < 0) continue;
            const int j = owner(bestCode);
            const int old = code[i];
            code[i] = bestCode;
            owners[bestCode] = i;
            if (j >= 0) {
                owners[old] = j;
                code[j] = old;
            } else {
                owners.erase(old);
            }
            improved = true;
        }
        if (!improved) break;
    }
    for (int i = 0; i < n; i++) codes[order[i]] = code[i];
    return true;
}

// A one-bit memory of the inputs and outputs: a flag that an input sets and
// an output resets (next = (flag or set) and not reset, or, if the set wins,
// set or (flag and not reset); set and reset are literals of an input and an
// output), or the previous value of an input or of an output.
struct Monitor {
    enum Kind { ResetWins, SetWins, PreviousInput, PreviousOutput };
    Kind kind;
    int input;
    bool inputValue;
    int output;
    bool outputValue;

    bool next(bool flag, const std::vector<uint8_t>& in, const std::vector<uint8_t>& out) const
    {
        switch (kind) {
        case PreviousInput: return in[input] != 0;
        case PreviousOutput: return out[output] != 0;
        default: break;
        }
        const bool set = (in[input] != 0) == inputValue;
        const bool reset = (out[output] != 0) == outputValue;
        return kind == SetWins ? set || (flag && !reset) : (flag || set) && !reset;
    }
};

// The monitors to try: flags for every input and output (and polarities),
// and the previous values of every input and output.
std::vector<Monitor> monitors(int uapCount, int capCount)
{
    std::vector<Monitor> result;
    for (int j = 0; j < capCount; j++) result.push_back(Monitor{Monitor::PreviousOutput, 0, true, j, true});
    for (int u = 0; u < uapCount; u++) result.push_back(Monitor{Monitor::PreviousInput, u, true, 0, true});
    for (int u = 0; u < uapCount; u++) {
        for (int j = 0; j < capCount; j++) {
            for (int variant = 0; variant < 8; variant++) {
                const Monitor::Kind kind = (variant & 4) != 0 ? Monitor::SetWins : Monitor::ResetWins;
                result.push_back(Monitor{kind, u, (variant & 1) != 0, j, (variant & 2) != 0});
            }
        }
    }
    return result;
}

// The product of the controller with the monitor (whose flag starts at 0):
// its states are the pairs (controller state, flag) that it reaches, the
// initial one first. Empty if the next flag is not determined by the letter.
Controller withMonitor(const Controller& ctl, const Monitor& mon,
                       const std::vector<std::vector<std::vector<uint8_t>>>& inputs, std::vector<int>& flagOf)
{
    std::map<std::pair<int, int>, int> index{{{ctl.initial, 0}, 0}};
    std::vector<std::pair<int, int>> pairs{{ctl.initial, 0}};
    Controller product;
    for (size_t i = 0; i < pairs.size(); i++) {
        const auto [c, f] = pairs[i];
        product.successor.emplace_back(inputs.size());
        product.output.push_back(ctl.output[c]);
        for (size_t l = 0; l < inputs.size(); l++) {
            int next = -1;
            for (const std::vector<uint8_t>& in : inputs[l]) {
                const int g = mon.next(f != 0, in, ctl.output[c][l]) ? 1 : 0;
                if (next != -1 && next != g) return Controller{};
                next = g;
            }
            const auto target = std::make_pair(ctl.successor[c][l], next);
            const auto inserted = index.emplace(target, static_cast<int>(pairs.size()));
            if (inserted.second) pairs.push_back(target);
            product.successor[i][l] = inserted.first->second;
        }
    }
    flagOf.clear();
    for (const auto& pair : pairs) flagOf.push_back(pair.second);
    return product;
}

// Codes for the states of a controller that contain the flags of monitors:
// such a bit only needs the input, the output and its own value to compute
// its next value. A monitor is taken if its flag tells states with the same
// flags apart and is determined by the controller state, possibly after
// splitting a few controller states by the flag (then ctl becomes the product
// with the monitor). The remaining bits number the states with the same
// flags, in breadth-first order (so the initial one gets code 0). Empty if
// there is no such monitor or too much to try.
std::vector<int> monitorCodes(Controller& ctl, const std::vector<std::vector<std::vector<uint8_t>>>& inputs,
                              int capCount, std::vector<int>& order, int& statebits, std::vector<Monitor>& chosen)
{
    const int uapCount = inputs.empty() || inputs[0].empty() ? 0 : static_cast<int>(inputs[0][0].size());
    long work = 0;
    for (const auto& letter : inputs) work += static_cast<long>(letter.size());
    const std::vector<Monitor> candidates = monitors(uapCount, capCount);
    work *= static_cast<long>(ctl.successor.size()) * static_cast<long>(candidates.size());
    if (work == 0 || work > MAX_MONITOR_WORK) return {};

    const int reachable = static_cast<int>(reachableStates(ctl).size());
    const int extra = std::max(1, reachable / MONITOR_SPLIT_SHARE);   // states that splitting may add, in total
    std::vector<int> cls(ctl.successor.size(), 0);    // the flags taken so far, as a number
    int classes = 1;
    int flags = 0;
    chosen.clear();
    std::vector<int> flagOf;
    for (const Monitor& mon : candidates) {
        if (flags >= 30) break;   // the codes are ints
        Controller product = withMonitor(ctl, mon, inputs, flagOf);
        if (product.successor.empty()) continue;
        const int current = static_cast<int>(reachableStates(ctl).size());
        const int added = static_cast<int>(product.successor.size()) - current;
        if (added > extra - (current - reachable)) continue;
        // the states of the product, with the flags taken so far (the controller state's) and this one
        std::vector<int> origin(product.successor.size());
        {
            std::map<std::pair<int, int>, int> index{{{ctl.initial, 0}, 0}};
            std::vector<std::pair<int, int>> pairs{{ctl.initial, 0}};
            for (size_t i = 0; i < pairs.size(); i++) {
                origin[i] = pairs[i].first;
                for (size_t l = 0; l < inputs.size(); l++) {
                    const auto target = std::make_pair(ctl.successor[pairs[i].first][l],
                                                       flagOf[product.successor[i][l]]);
                    if (index.emplace(target, static_cast<int>(pairs.size())).second) pairs.push_back(target);
                }
            }
        }
        std::set<std::pair<int, int>> split;
        for (size_t i = 0; i < origin.size(); i++) split.emplace(cls[origin[i]], flagOf[i]);
        if (static_cast<int>(split.size()) == classes) continue;
        std::vector<int> next(product.successor.size());
        for (size_t i = 0; i < origin.size(); i++) next[i] = (cls[origin[i]] << 1) | flagOf[i];
        cls = std::move(next);
        classes = static_cast<int>(split.size());
        ctl = std::move(product);
        chosen.push_back(mon);
        flags++;
    }
    if (flags == 0) return {};
    order = reachableStates(ctl);
    std::map<int, int> used;
    std::vector<int> index(ctl.successor.size(), 0);
    int largest = 1;
    for (int c : order) {
        index[c] = used[cls[c]]++;
        largest = std::max(largest, used[cls[c]]);
    }
    int rest = 0;
    while ((1 << rest) < largest) rest++;
    statebits = std::max(1, flags + rest);
    if (statebits > 30) return {};
    std::vector<int> codes(ctl.successor.size(), -1);
    for (int c : order) codes[c] = (cls[c] << rest) | index[c];
    return codes;
}

// Restrict the moves to outputs that keep the flags in the codes (made by
// monitorCodes, with the given number of further bits) monitor flags: from
// every state, for every input of the letter, the monitor computes the flag
// of the successor. The codes are only labels, so this only makes the
// circuit smaller; the deterministic moves of the controller always remain.
Moves keepFlags(const Moves& moves, const std::vector<Monitor>& chosen, const std::vector<int>& codes, int rest,
                const std::vector<std::vector<std::vector<uint8_t>>>& inputs, int capCount, OutputTable& table,
                MTBDD capVars)
{
    const int flags = static_cast<int>(chosen.size());
    Moves result(moves.size());
    MTBDD set = mtbdd_false, var = mtbdd_false;
    mtbdd_refs_pushptr(&set);
    mtbdd_refs_pushptr(&var);
    std::vector<uint8_t> out(capCount, 0);
    for (size_t c = 0; c < moves.size(); c++) {
        if (codes[c] == -1 || moves[c].empty()) continue;
        result[c].resize(moves[c].size());
        for (size_t l = 0; l < moves[c].size(); l++) {
            for (const Choice& choice : moves[c][l]) {
                set = table.sets()[choice.outputs];
                bool ok = true;
                for (int t = 0; t < flags && ok; t++) {
                    const Monitor& mon = chosen[t];
                    const int bit = rest + flags - 1 - t;
                    const bool flag = ((codes[c] >> bit) & 1) != 0;
                    const bool target = ((codes[choice.successor] >> bit) & 1) != 0;
                    bool can[2];
                    for (int v = 0; v < 2; v++) {
                        if (mon.kind != Monitor::PreviousInput) out[mon.output] = static_cast<uint8_t>(v);
                        can[v] = true;
                        for (const std::vector<uint8_t>& in : inputs[l]) can[v] = can[v] && mon.next(flag, in, out) == target;
                    }
                    if (mon.kind != Monitor::PreviousInput) out[mon.output] = 0;
                    if (can[0] && can[1]) continue;
                    if (!can[0] && !can[1]) {
                        ok = false;
                        break;
                    }
                    var = sylvan_ithvar(mtbdd_set_first(capVars) + static_cast<uint32_t>(mon.output));
                    set = sylvan_and(set, can[1] ? var : sylvan_not(var));
                    ok = set != mtbdd_false;
                }
                if (ok) result[c][l].push_back(Choice{table.add(set), choice.successor});
            }
        }
    }
    mtbdd_refs_popptr(2);
    return result;
}

// The cube of a state code over statevars (most significant bit first)
// without the bits in free, then priority 0 over priovars (if any).
MTBDD codeCube(uint32_t code, uint32_t free, MTBDD statevars, MTBDD priovars)
{
    std::vector<uint32_t> vars;
    for (MTBDD v = statevars; v != mtbdd_set_empty(); v = mtbdd_set_next(v)) vars.push_back(mtbdd_set_first(v));
    MTBDD cube = mtbdd_true;
    mtbdd_refs_pushptr(&cube);
    for (size_t bit = 0; bit < vars.size(); bit++) {
        const uint32_t var = vars[vars.size() - 1 - bit];
        if ((free >> bit) & 1) continue;
        cube = ((code >> bit) & 1) ? mtbdd_makenode(var, mtbdd_false, cube) : mtbdd_makenode(var, cube, mtbdd_false);
    }
    vars.clear();
    for (MTBDD v = priovars; v != mtbdd_set_empty(); v = mtbdd_set_next(v)) vars.push_back(mtbdd_set_first(v));
    while (!vars.empty()) {
        cube = mtbdd_makenode(vars.back(), cube, mtbdd_false);
        vars.pop_back();
    }
    mtbdd_refs_popptr(1);
    return cube;
}

// The game whose states are the controller states reachable from the initial
// one (in the given order, the initial one first, with the given codes; the
// initial one has code 0), with exactly the controller's moves. With free (by
// state), a state's code is a cube: the bits in its free mask take any value,
// so its successors may set them freely, and every code of the cube behaves
// as the state.
std::unique_ptr<SymGame> toGame(SymGame& g, const std::vector<int>& order, const std::vector<int>& number,
                                int statebits, const Moves& moves, const std::vector<MTBDD>& letterSets,
                                const std::vector<MTBDD>& sets, const std::vector<uint32_t>& free = {})
{
    auto freeOf = [&](int state) { return free.empty() ? 0u : free[state]; };
    const int count = static_cast<int>(order.size());
    auto res = std::make_unique<SymGame>(statebits, g.priobits, g.uap_count, g.cap_count, g.maxprio);

    BddPool pool;
    // rename the uap and cap variables of the letters and outputs to the new game
    MTBDD rename = pool.keep(mtbdd_map_empty());
    auto renameVars = [&](MTBDD fromVars, MTBDD toVars, int vars) {
        for (int j = 0; j < vars; j++) {
            const uint32_t from = mtbdd_set_first(fromVars) + static_cast<uint32_t>(j);
            const uint32_t to = mtbdd_set_first(toVars) + static_cast<uint32_t>(j);
            const MTBDD var = pool.keep(sylvan_ithvar(to));
            rename = pool.keep(mtbdd_map_add(rename, from, var));
        }
    };
    renameVars(g.uap_vars, res->uap_vars, g.uap_count);
    renameVars(g.cap_vars, res->cap_vars, g.cap_count);
    std::vector<MTBDD> letters;
    for (MTBDD set : letterSets) letters.push_back(pool.keep(sylvan_compose(set, rename)));
    std::map<int, MTBDD> outputs;

    MTBDD state = mtbdd_false, input = mtbdd_false, move = mtbdd_false, target = mtbdd_false, step = mtbdd_false;
    for (MTBDD* r : {&state, &input, &move, &target, &step}) mtbdd_refs_pushptr(r);
    for (int i = 0; i < count; i++) {
        state = codeCube(static_cast<uint32_t>(number[order[i]]), freeOf(order[i]), res->s_vars, mtbdd_set_empty());
        for (size_t l = 0; l < letters.size(); l++) {
            for (const Choice& choice : moves[order[i]][l]) {
                auto out = outputs.find(choice.outputs);
                if (out == outputs.end()) {
                    out = outputs.emplace(choice.outputs, pool.keep(sylvan_compose(sets[choice.outputs], rename))).first;
                }
                input = sylvan_and(state, letters[l]);
                move = sylvan_and(input, out->second);
                target = codeCube(static_cast<uint32_t>(number[choice.successor]), freeOf(choice.successor), res->ns_vars,
                                  res->np_vars);
                step = sylvan_and(move, target);
                res->strategies = sylvan_or(res->strategies, move);
                res->trans = sylvan_or(res->trans, step);
            }
        }
    }
    mtbdd_refs_popptr(5);
    return res;
}

// Fill in the result for a correct controller of the machine of ex: the game
// with exactly the controller's moves and, with relax, variants for the
// encoder: with all moves that keep it correct (in breadth-first and
// depth-first numbering), and with codes of monitor flags for the product of
// the controller with the monitors (with its moves, and with the moves that
// keep it correct and keep the flags). Unless the moves of the machine are a
// winning multi-strategy (with checkCycles), relaxed moves are only used if
// no play under them is lost. False if the controller makes a move the
// machine does not allow.
// A state bit that is an output bit: of the moves into the state ("incoming": the
// next-state function of the latch is then the function of the output, which the
// circuit already has) or of the moves out of it ("outgoing": the output is then
// the latch itself), stored inverted if needed so that the initial state has code 0.
struct OutputBit {
    int output;
    bool incoming;
    bool invert;
};

// State codes whose leading bits are output bits, chosen with all correct moves:
// a state can take such a bit if every move into it (or out of it) can give the
// output that value. Bits are taken greedily while they tell apart states with
// the same bits so far; the remaining bits number the states with the same
// leading bits (breadth-first). Returns the codes (by state, -1 if unreachable),
// the order of the states, the number of bits and the output bits, or no codes
// if no output bit tells states apart.
std::vector<int> outputCodes(const Moves& all, int initial, int capCount, const std::vector<MTBDD>& sets,
                             MTBDD capVars, std::vector<int>& order, int& statebits, std::vector<OutputBit>& chosen)
{
    order = numbering(initial, all, false);
    const size_t n = all.size();
    std::vector<char> reachable(n, 0);
    for (int c : order) reachable[c] = 1;
    std::vector<int> cls(n, 0);   // the leading bits so far
    int classes = 1;
    chosen.clear();
    MTBDD var = mtbdd_false, test = mtbdd_false;
    mtbdd_refs_pushptr(&var);
    mtbdd_refs_pushptr(&test);
    for (int j = 0; j < capCount && chosen.size() < 30; j++) {
        var = sylvan_ithvar(mtbdd_set_first(capVars) + static_cast<uint32_t>(j));
        // per state, the values of output j that every move into it (out of it) allows: bit 0 for
        // the value 0, bit 1 for the value 1
        std::vector<int> into(n, 3), outOf(n, 3);
        for (size_t c = 0; c < n; c++) {
            if (!reachable[c]) continue;
            for (const std::vector<Choice>& letter : all[c]) {
                for (const Choice& choice : letter) {
                    int values = 0;
                    test = sylvan_and(sets[choice.outputs], sylvan_not(var));
                    if (test != mtbdd_false) values |= 1;
                    test = sylvan_and(sets[choice.outputs], var);
                    if (test != mtbdd_false) values |= 2;
                    into[choice.successor] &= values;
                    outOf[c] &= values;
                }
            }
        }
        for (const bool incoming : {false, true}) {
            if (chosen.size() >= 30) break;
            const std::vector<int>& allowed = incoming ? into : outOf;
            bool possible = true;
            for (int c : order) possible = possible && allowed[c] != 0;
            if (!possible) continue;
            const bool invert = allowed[initial] == 2;
            // the forced values first; each free state takes the value that tells the most apart
            std::vector<int> value(n, 0);
            for (int c : order) value[c] = allowed[c] == 3 ? 0 : (allowed[c] == 2 ? 1 : 0) ^ (invert ? 1 : 0);
            auto count = [&]() {
                std::set<std::pair<int, int>> split;
                for (int c : order) split.emplace(cls[c], value[c]);
                return static_cast<int>(split.size());
            };
            int best = count();
            for (int c : order) {
                if (c == initial || allowed[c] != 3) continue;
                value[c] = 1;
                const int with = count();
                if (with > best) {
                    best = with;
                } else {
                    value[c] = 0;
                }
            }
            if (best == classes) continue;
            for (int c : order) cls[c] = (cls[c] << 1) | value[c];
            classes = best;
            chosen.push_back(OutputBit{j, incoming, invert});
        }
    }
    mtbdd_refs_popptr(2);
    if (chosen.empty()) return {};
    const int flags = static_cast<int>(chosen.size());
    std::map<int, int> used;
    std::vector<int> index(n, 0);
    int largest = 1;
    for (int c : order) {
        index[c] = used[cls[c]]++;
        largest = std::max(largest, used[cls[c]]);
    }
    int rest = 0;
    while ((1 << rest) < largest) rest++;
    statebits = std::max(1, flags + rest);
    if (statebits > 30) return {};   // the codes are ints
    std::vector<int> codes(n, -1);
    for (int c : order) codes[c] = (cls[c] << rest) | index[c];
    return codes;
}

// Restrict the moves to the outputs that agree with the output bits of the codes.
Moves keepOutputBits(const Moves& moves, const std::vector<OutputBit>& chosen, const std::vector<int>& codes, int rest,
                     OutputTable& table, MTBDD capVars)
{
    const int flags = static_cast<int>(chosen.size());
    Moves result(moves.size());
    MTBDD set = mtbdd_false, var = mtbdd_false;
    mtbdd_refs_pushptr(&set);
    mtbdd_refs_pushptr(&var);
    for (size_t c = 0; c < moves.size(); c++) {
        if (codes[c] == -1 || moves[c].empty()) continue;
        result[c].resize(moves[c].size());
        for (size_t l = 0; l < moves[c].size(); l++) {
            for (const Choice& choice : moves[c][l]) {
                set = table.sets()[choice.outputs];
                bool possible = true;
                for (int t = 0; t < flags && possible; t++) {
                    const OutputBit& b = chosen[t];
                    const int bit = rest + flags - 1 - t;
                    const int stored = (codes[b.incoming ? choice.successor : static_cast<int>(c)] >> bit) & 1;
                    const bool value = (stored ^ (b.invert ? 1 : 0)) != 0;
                    var = sylvan_ithvar(mtbdd_set_first(capVars) + static_cast<uint32_t>(b.output));
                    set = sylvan_and(set, value ? var : sylvan_not(var));
                    possible = set != mtbdd_false;
                }
                if (possible) result[c][l].push_back(Choice{table.add(set), choice.successor});
            }
        }
    }
    mtbdd_refs_popptr(2);
    return result;
}

bool emit(SymGame& game, const Extraction& ex, const Controller& ctl, const MealyMinimizer::Options& options,
          bool checkCycles, BddPool& pool, MealyMinimizer::Result& result)
{
    const MealyMinimizer::Variants& kinds = options.variants;
    OutputTable table(pool);
    const Moves moves = movesOf(ctl, table, game.cap_vars);
    const std::vector<std::vector<char>> pairs = reached(ex.machine, ex.outputSets, table.sets(), ctl.initial, moves);
    if (pairs.empty()) return false;
    const std::vector<int> order = numbering(ctl.initial, moves, false);
    int statebits;
    std::vector<int> codes = positions(order, moves.size(), statebits);
    result.game = toGame(game, order, codes, statebits, moves, ex.letterSets, table.sets());
    result.states = static_cast<int>(order.size());
    result.variants.clear();
    if (!options.relax) return true;

    auto correct = [&](int initial, const Moves& all) {
        return !reached(ex.machine, ex.outputSets, table.sets(), initial, all).empty()
               && (!checkCycles || winningMoves(ex.machine, ex.outputSets, table.sets(), initial, all));
    };
    const Moves all = kinds.relaxed ? relaxedMoves(ex.machine, ex.outputSets, moves.size(), pairs, table) : moves;
    if (correct(ctl.initial, all)) {
        if (kinds.relaxed && kinds.renumber) {
            // the numbering of the states matters for the encoding, in no predictable way
            for (const bool depthFirst : {false, true}) {
                const std::vector<int> numbers = numbering(ctl.initial, all, depthFirst);
                codes = positions(numbers, all.size(), statebits);
                result.variants.push_back({toGame(game, numbers, codes, statebits, all, ex.letterSets, table.sets()),
                                           depthFirst ? "relaxed, depth-first" : "relaxed, breadth-first",
                                           MealyMinimizer::Kind::relaxed});
            }
        } else if (kinds.relaxed) {
            codes = positions(order, all.size(), statebits);
            result.variants.push_back({toGame(game, order, codes, statebits, all, ex.letterSets, table.sets()),
                                       "relaxed", MealyMinimizer::Kind::relaxed});
        }
        if (kinds.relaxed && kinds.cubeCodes) {
            for (const bool byIncoming : {false, true}) {
                std::vector<uint32_t> free;
                const std::vector<int> numbers = cubeNumbering(ctl.initial, all, byIncoming, codes, statebits, free);
                if (numbers.empty()) break;
                result.variants.push_back(
                    {toGame(game, numbers, codes, statebits, all, ex.letterSets, table.sets(), free),
                     byIncoming ? "relaxed, cube codes by incoming moves" : "relaxed, cube codes", MealyMinimizer::Kind::cubeCodes});
            }
        }
        if (kinds.relaxed && kinds.adjacencyCodes) {
            const std::vector<int> numbers = numbering(ctl.initial, all, false);
            codes = positions(numbers, all.size(), statebits);
            if (improveCodes(all, numbers, codes, statebits, statebits)) {
                result.variants.push_back(
                    {toGame(game, numbers, codes, statebits, all, ex.letterSets, table.sets()), "relaxed, adjacency codes",
                     MealyMinimizer::Kind::adjacencyCodes});
            }
        }
        // state bits that are output bits cost (almost) no gates
        std::vector<int> outputOrder;
        std::vector<OutputBit> outputBits;
        codes.clear();
        if (kinds.outputCodes) {
            codes = outputCodes(all, ctl.initial, game.cap_count, table.sets(), game.cap_vars, outputOrder, statebits,
                                outputBits);
        }
        if (!codes.empty()) {
            const int flags = static_cast<int>(outputBits.size());
            const Moves agreeing = keepOutputBits(all, outputBits, codes, statebits - flags, table, game.cap_vars);
            if (correct(ctl.initial, agreeing)) {
                // only the states that these moves reach (where correct() checked that every state has a
                // move for every letter); the codes are labels, and the initial state keeps code 0
                const std::vector<int> agreeingOrder = numbering(ctl.initial, agreeing, false);
                if (options.verbose) {
                    std::cerr << "state codes of " << flags << " output bits and " << statebits - flags
                              << " further bits for " << agreeingOrder.size() << " states" << std::endl;
                }
                result.variants.push_back(
                    {toGame(game, agreeingOrder, codes, statebits, agreeing, ex.letterSets, table.sets()), "output bits", MealyMinimizer::Kind::outputCodes});
                // with the further bits by adjacency
                if (kinds.adjacencyCodes && improveCodes(agreeing, agreeingOrder, codes, statebits, statebits - flags)) {
                    result.variants.push_back({toGame(game, agreeingOrder, codes, statebits, agreeing, ex.letterSets,
                                                      table.sets()), "output bits, adjacency codes", MealyMinimizer::Kind::outputCodes});
                }
            }
        }
    }

    if (!kinds.monitorCodes) return true;
    Controller split = ctl;
    std::vector<int> splitOrder;
    std::vector<Monitor> chosen;
    const auto inputs = ex.inputs();
    codes = monitorCodes(split, inputs, game.cap_count, splitOrder, statebits, chosen);
    if (codes.empty()) return true;
    const Moves exact = movesOf(split, table, game.cap_vars);
    const std::vector<std::vector<char>> splitPairs =
        reached(ex.machine, ex.outputSets, table.sets(), split.initial, exact);
    if (splitPairs.empty()) return true;
    const int flags = static_cast<int>(chosen.size());
    if (options.verbose) {
        std::cerr << "state codes of " << flags << " monitor flags and " << statebits - flags << " further bits for "
                  << splitOrder.size() << " states" << std::endl;
    }
    // with the controller's moves, and with the freedom that keeps the flags
    result.variants.push_back({toGame(game, splitOrder, codes, statebits, exact, ex.letterSets, table.sets()),
                               "monitor codes", MealyMinimizer::Kind::monitorCodes});
    const size_t exactVariant = result.variants.size() - 1;
    if (!kinds.relaxed) return true;
    const Moves free = keepFlags(relaxedMoves(ex.machine, ex.outputSets, exact.size(), splitPairs, table), chosen,
                                 codes, statebits - flags, inputs, game.cap_count, table, game.cap_vars);
    if (correct(split.initial, free)) {
        result.variants[exactVariant].covered = true;
        result.variants.push_back({toGame(game, splitOrder, codes, statebits, free, ex.letterSets, table.sets()),
                                   "monitor codes, relaxed", MealyMinimizer::Kind::relaxedMonitorCodes});
        // with the further bits by adjacency
        if (kinds.adjacencyCodes && improveCodes(free, splitOrder, codes, statebits, statebits - flags)) {
            result.variants.push_back({toGame(game, splitOrder, codes, statebits, free, ex.letterSets, table.sets()),
                                       "monitor codes, relaxed, adjacency codes", MealyMinimizer::Kind::relaxedMonitorCodes});
        }
    }
    return true;
}

// Search for a small controller that may make any of the given moves, a
// superset of the moves of the strategy. If they are a winning multi-strategy
// (such as the permissive extension), the search starts from the controller of
// their bisimulation quotient, as minimize() does; otherwise (such as all moves
// into the winning region) it only looks for controllers with fewer than bound
// states, and every candidate of the SAT solver is checked for losing cycles,
// which are then excluded by rank constraints. The lower bound of these moves
// is weaker than that of the strategy (more states are compatible), so the
// encoding grows with the square of the number of states: search upwards from
// the lower bound (by one state when the candidates are checked for losing
// cycles, see StateSearch), and bisect once a controller is found.
// Fills in the result if a controller with fewer than bound states is found.
// With relax, a controller that is not smaller still gives games with its
// freedom to more, for the encoder to choose from (fewer states need not give
// a smaller circuit).
bool searchWider(SymGame& game, MTBDD moves, const char* name, int bound, long work,
                 const MealyMinimizer::Options& options, MealyMinimizer::Result& result,
                 std::vector<MealyMinimizer::Variant>& more, std::vector<MealyMinimizer::Search>& searches)
{
    const bool verbose = options.verbose;
    const bool relax = options.relax;
    const size_t index = searches.size();
    searches.emplace_back();
    searches[index].moves = name;
    BddPool pool;
    Extraction ex(game, moves, pool);
    try {
        ex.run();
    } catch (const Aborted& a) {
        if (verbose) std::cerr << "search over the " << name << " skipped: " << a.reason << std::endl;
        searches[index].skipped = a.reason;
        return false;
    }
    const bool checkCycles = !winning(ex.machine);
    std::vector<int> block;
    const Machine q = quotient(ex.machine, block, checkCycles);
    Overlap overlap(ex.outputSets);
    std::vector<std::vector<char>> compatible;
    const std::vector<int> anchors = incompatibleStates(q, overlap, options.strongBound,
                                                        options.winningBound && checkCycles, compatible);
    const int lowerBound = static_cast<int>(anchors.size());
    searches[index].machine = ex.machine.states();
    searches[index].quotient = q.states();
    searches[index].lowerBound = lowerBound;
    if (verbose) {
        std::cerr << name << ": " << ex.machine.states() << " states, " << q.states()
                  << " after bisimulation, at least " << lowerBound << " states needed" << std::endl;
    }

    const uint32_t firstCap = game.cap_count > 0 ? mtbdd_set_first(game.cap_vars) : 0;
    Controller best;
    bool haveController = false;            // best is a controller of these moves
    int known = bound;                      // the number of states below which to search
    if (!checkCycles) {
        best = fromQuotient(q, ex.outputSets, game);
        known = static_cast<int>(reachableStates(best).size());
        haveController = true;
    }
    StateSearch sizes(lowerBound, known, options.fewerFirst, checkCycles && options.unitSteps);
    long remaining = work;
    for (int k = sizes.next(); k > 0 && remaining > 0; k = sizes.next()) {
        Encoding enc(q, ex.outputSets, game.cap_count, firstCap, anchors, compatible, k, checkCycles);
        if (enc.estimatedSize() > MAX_CLAUSES) {
            if (verbose) std::cerr << name << ": the encoding for " << k << " states is too large" << std::endl;
            searches[index].calls.push_back({k, -2, 0});
            sizes.tooLarge(k);
            continue;
        }
        // every attempt may use half of the remaining work, with its clauses and the rounds of
        // candidates and exclusions
        const long allowance = (remaining + 1) / 2;
        int r;
        long cycles = 0;
        Controller candidate;
        try {
            r = enc.build(allowance);
            if (r == l_True) r = enc.solve(allowance);
            while (r == l_True) {
                candidate = enc.controller();
                if (!checkCycles) break;
                const std::vector<Step> cycle = losingCycle(ex.machine, block, ex.outputSets, candidate, firstCap);
                if (cycle.empty()) break;
                cycles++;
                r = enc.exclude(cycle) ? enc.solve(allowance) : l_False;
            }
        } catch (const Aborted&) {
            r = l_Undef;
        }
        remaining -= enc.work();
        if (verbose) {
            // every correct controller only makes moves into the winning region, so over all
            // winning moves, "do not suffice" is a lower bound for every controller of the specification
            std::cerr << name << ": " << k << " states "
                      << (r == l_True ? "suffice" : r == l_False ? "do not suffice" : "undecided") << " ("
                      << enc.variables() << " variables, " << enc.clauses() << " clauses, work " << enc.work()
                      << ", " << cycles << " losing cycles excluded)" << std::endl;
        }
        searches[index].calls.push_back({k, r == l_True ? 1 : r == l_False ? 0 : -1, enc.work()});
        if (r == l_True) {
            best = std::move(candidate);
            sizes.suffice(static_cast<int>(reachableStates(best).size()));
            haveController = true;
        } else if (r == l_False) {
            sizes.refute(k);
        } else {
            sizes.undecide(k);
        }
    }
    searches[index].stop = sizes.proven() ? "proven" : remaining <= 0 ? "work" : "undecided";
    searches[index].proven = sizes.proven();
    if (!haveController) return false;
    const int top = sizes.states();
    searches[index].found = top;
    if (top >= bound) {
        // not smaller than the controller of the strategy, but with more freedom
        MealyMinimizer::Result other;
        if (relax && !checkCycles && emit(game, ex, best, options, false, pool, other)) {
            const std::string prefix = std::string(name) + ": ";
            more.push_back({std::move(other.game), prefix + "controller", MealyMinimizer::Kind::controller, false});
            for (auto& variant : other.variants) {
                more.push_back({std::move(variant.game), prefix + variant.name, variant.kind, false, variant.covered});
            }
        }
        return false;
    }

    // check the controller independently of the encoding before using it
    bool correct;
    try {
        correct = losingCycle(ex.machine, block, ex.outputSets, best, firstCap).empty();
    } catch (const Aborted&) {
        correct = false;
    }
    MealyMinimizer::Result found;
    if (!correct || !emit(game, ex, best, options, checkCycles, pool, found)) {
        // cannot happen; keep the other controller rather than risk a wrong one
        std::cerr << "error: the controller found over the " << name << " is not winning" << std::endl;
        return false;
    }
    result = std::move(found);
    if (verbose) std::cerr << "a controller with " << result.states << " states uses the " << name << std::endl;
    return true;
}

}  // namespace

const char* MealyMinimizer::name(Kind kind)
{
    switch (kind) {
        case Kind::controller: return "controller";
        case Kind::relaxed: return "relaxed";
        case Kind::cubeCodes: return "cube codes";
        case Kind::adjacencyCodes: return "adjacency codes";
        case Kind::outputCodes: return "output codes";
        case Kind::monitorCodes: return "monitor codes";
        case Kind::relaxedMonitorCodes: return "relaxed monitor codes";
        case Kind::unclarified: return "unclarified";
    }
    return "";
}

MealyMinimizer::Result MealyMinimizer::minimize(SymGame& game, const Options& options)
{
    const bool verbose = options.verbose;
    const bool relax = options.relax;
    Result result;
    std::vector<Search> searches(1);
    searches.reserve(3);   // the strategy, the extended strategy and all winning moves: search stays valid
    Search& search = searches[0];
    search.moves = "strategy";
    BddPool pool;
    Extraction ex(game, game.strategies, pool);
    try {
        ex.run();
    } catch (const Aborted& a) {
        if (verbose) std::cerr << "mealy minimisation skipped: " << a.reason << std::endl;
        search.skipped = a.reason;
        result.searches = std::move(searches);
        return result;
    }
    if (!winning(ex.machine)) {
        // the freedom in the strategy is not a winning multi-strategy (e.g. the moves that
        // the symbolic solver keeps); minimising over it could give a wrong controller
        if (verbose) std::cerr << "mealy minimisation skipped: the strategy allows losing cycles" << std::endl;
        search.skipped = "the strategy allows losing cycles";
        result.searches = std::move(searches);
        return result;
    }

    std::vector<int> block;
    const Machine q = quotient(ex.machine, block, false);
    Overlap overlap(ex.outputSets);
    std::vector<std::vector<char>> compatible;
    const std::vector<int> anchors = incompatibleStates(q, overlap, options.strongBound, false, compatible);
    const int lowerBound = static_cast<int>(anchors.size());
    search.machine = ex.machine.states();
    search.quotient = q.states();
    search.lowerBound = lowerBound;
    if (verbose) {
        std::cerr << "mealy machine: " << ex.machine.states() << " states, " << q.states() << " after bisimulation, "
                  << ex.letterSets.size() << " input classes, at least " << lowerBound << " states needed"
                  << std::endl;
    }

    // The bisimulation quotient is a controller; look for smaller ones, first
    // with as many states as the lower bound (which often suffices), as
    // StateSearch chooses. Every attempt may use half of the remaining work.
    const uint32_t firstCap = game.cap_count > 0 ? mtbdd_set_first(game.cap_vars) : 0;
    Controller best = fromQuotient(q, ex.outputSets, game);
    StateSearch sizes(lowerBound, static_cast<int>(reachableStates(best).size()), options.fewerFirst);
    long remaining = options.work;
    for (int k = sizes.next(); k > 0 && remaining > 0; k = sizes.next()) {
        Encoding enc(q, ex.outputSets, game.cap_count, firstCap, anchors, compatible, k, false);
        if (enc.estimatedSize() > MAX_CLAUSES) {
            if (verbose) std::cerr << "mealy minimisation: the encoding for " << k << " states is too large" << std::endl;
            search.calls.push_back({k, -2, 0});
            sizes.tooLarge(k);
            continue;
        }
        const long allowance = (remaining + 1) / 2;
        int r = enc.build(allowance);
        if (r == l_True) r = enc.solve(allowance);
        remaining -= enc.work();
        if (verbose) {
            std::cerr << "mealy minimisation: " << k << " states "
                      << (r == l_True ? "suffice" : r == l_False ? "do not suffice" : "undecided") << " ("
                      << enc.variables() << " variables, " << enc.clauses() << " clauses, work " << enc.work()
                      << ")" << std::endl;
        }
        search.calls.push_back({k, r == l_True ? 1 : r == l_False ? 0 : -1, enc.work()});
        if (r == l_True) {
            best = enc.controller();
            sizes.suffice(static_cast<int>(reachableStates(best).size()));
        } else if (r == l_False) {
            sizes.refute(k);
        } else {
            sizes.undecide(k);
        }
    }
    if (verbose && !sizes.proven()) {
        std::cerr << "mealy minimisation: " << sizes.states() << " states, not proven minimal for this strategy" << std::endl;
    }
    search.found = sizes.states();
    search.proven = sizes.proven();
    search.stop = sizes.proven() ? "proven" : remaining <= 0 ? "work" : "undecided";
    if (!emit(game, ex, best, options, false, pool, result)) {
        // cannot happen; keep the original controller rather than risk a wrong one
        std::cerr << "error: the minimised controller makes moves the strategy does not allow" << std::endl;
        Result none;
        none.searches = std::move(searches);
        return none;
    }

    // then look for smaller controllers with more freedom; with relax, the games of the
    // controllers that are not used are variants, as fewer states need not give a smaller circuit
    std::vector<Variant> more;
    std::string current = "strategy";   // the moves of the controller in result
    struct Wider {
        MTBDD moves;
        const char* name;
        long work;
        int maxStates;   // only search if the controller has at most this many states (0: always)
    };
    for (const Wider& w : {Wider{options.extended, "extended strategy", options.work, 0},
                           Wider{options.anyMoves, "all winning moves", options.anyWork, options.anyStates}}) {
        Result smaller;
        if (w.moves == mtbdd_false) continue;
        if (w.maxStates > 0 && result.states > w.maxStates) {
            searches.emplace_back();
            searches.back().moves = w.name;
            searches.back().skipped = "the controller has more than " + std::to_string(w.maxStates) + " states";
            continue;
        }
        if (!searchWider(game, w.moves, w.name, result.states, w.work, options, smaller, more, searches)) continue;
        if (relax) {
            more.push_back({std::move(result.game), current + ": controller", Kind::controller, false});
            for (auto& variant : result.variants) {
                more.push_back({std::move(variant.game), current + ": " + variant.name, variant.kind, false,
                                variant.covered});
            }
        }
        result = std::move(smaller);
        current = w.name;
    }
    for (auto& variant : more) result.variants.push_back(std::move(variant));
    result.searches = std::move(searches);
    return result;
}
