/**
 * Copyright Tom van Dijk
 */

/*
 * Conversion of acceptance conditions to parity (paritizeHoa, see simplehoa.h).
 *
 * Whether a run of an automaton is accepting depends only on the set of
 * "atoms" that it sees infinitely often, where an atom is a combination of
 * acceptance sets that occurs on edges. The Zielonka tree of the acceptance
 * condition over the atoms has the set of all atoms as its root, and as the
 * children of a node the maximal subsets of the node with the other
 * acceptance. If the tree is a chain, every edge simply gets the priority of
 * the deepest node that contains its atom. Otherwise the automaton is
 * multiplied with the Zielonka tree automaton, whose states are the leaves of
 * the tree, which is a deterministic parity automaton for the condition with
 * the minimal number of states (Casares, Colcombet and Fijalkow, "Optimal
 * transformations of games and automata using Muller conditions", 2021).
 */

#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "hoainternal.hpp"
#include "simplehoa.h"

using namespace hoa;

namespace {

struct ConversionError {
    int code;
    std::string message;
};

[[noreturn]] void conversionError(std::string message, int code = 100)
{
    throw ConversionError{code, std::move(message)};
}

// The acceptance condition as a function of the set of atoms that a run sees
// infinitely often: Inf(x) holds iff an atom that contains x is seen, and
// Inf(!x) iff an atom without x is seen. Sets of acceptance sets and sets of
// atoms are bit masks.
class AtomCondition final {
public:
    // Translate acc over the given atoms; false if acc is not an acceptance
    // condition of the HOA format.
    bool build(const BTree* acc, const std::vector<uint64_t>& atoms)
    {
        nodes.clear();
        relevant = 0;
        root = translate(acc, atoms);
        return root >= 0;
    }

    bool accepts(uint64_t seen) const { return eval(root, 0, seen) == 1; }

    // Find a set of atoms X with the given acceptance, seen <= X <= seen |
    // open, that is not a subset of any mask in avoid. The witness is the
    // largest such X that the search ends on.
    bool find(uint64_t open, uint64_t seen, bool target, const std::vector<uint64_t>& avoid, uint64_t& witness)
    {
        for (uint64_t m : avoid) {
            if (((seen | open) & ~m) == 0) return false;
        }
        const int r = eval(root, open, seen);
        if (r != 2) {
            witness = seen | open;
            return r == (target ? 1 : 0);
        }
        if (--budget < 0) conversionError("the acceptance condition is too complex to analyse");
        const uint64_t branch = open & relevant;
        const uint64_t atom = branch & (~branch + 1);
        return find(open & ~atom, seen | atom, target, avoid, witness)
               || find(open & ~atom, seen, target, avoid, witness);
    }

private:
    struct Node {
        NodeType type;          // NT_BOOL, NT_AND, NT_OR, NT_INF or NT_FIN
        int left;
        int right;
        uint64_t atoms;         // NT_INF/NT_FIN: the atoms that satisfy the Inf
        bool value;             // NT_BOOL
    };

    // Bounds the recursion depth of translate and eval.
    static constexpr size_t MAX_NODES = 4096;

    int add(const Node& n)
    {
        if (nodes.size() == MAX_NODES) conversionError("the acceptance condition is too large");
        nodes.push_back(n);
        return static_cast<int>(nodes.size()) - 1;
    }

    int translate(const BTree* t, const std::vector<uint64_t>& atoms)
    {
        if (t == nullptr) return -1;
        switch (t->type) {
            case NT_BOOL:
                return add(Node{NT_BOOL, -1, -1, 0, t->id != 0});
            case NT_AND:
            case NT_OR: {
                const int l = translate(t->left, atoms);
                const int r = translate(t->right, atoms);
                if (l < 0 || r < 0) return -1;
                return add(Node{t->type, l, r, 0, false});
            }
            case NT_FIN:
            case NT_INF: {
                const BTree* set = t->left;
                const bool negated = set != nullptr && set->type == NT_NOT;
                if (negated) set = set->left;
                if (set == nullptr || set->type != NT_SET || set->id < 0 || set->id >= 64) return -1;
                uint64_t mask = 0;
                for (size_t i = 0; i < atoms.size(); i++) {
                    const bool member = (atoms[i] >> set->id) & 1;
                    if (member != negated) mask |= uint64_t(1) << i;
                }
                relevant |= mask;
                return add(Node{t->type, -1, -1, mask, false});
            }
            default:
                return -1;
        }
    }

    // Three-valued evaluation (0 false, 1 true, 2 unknown) for the sets of
    // atoms X with seen <= X <= seen | open.
    int eval(int n, uint64_t open, uint64_t seen) const
    {
        const Node& node = nodes[n];
        switch (node.type) {
            case NT_BOOL:
                return node.value ? 1 : 0;
            case NT_INF:
            case NT_FIN: {
                const int inf = (node.atoms & seen) != 0 ? 1 : (node.atoms & open) != 0 ? 2 : 0;
                return node.type == NT_INF || inf == 2 ? inf : 1 - inf;
            }
            case NT_AND: {
                const int l = eval(node.left, open, seen);
                if (l == 0) return 0;
                const int r = eval(node.right, open, seen);
                return r == 0 ? 0 : (l == 1 && r == 1 ? 1 : 2);
            }
            default: {
                const int l = eval(node.left, open, seen);
                if (l == 1) return 1;
                const int r = eval(node.right, open, seen);
                return r == 1 ? 1 : (l == 0 && r == 0 ? 0 : 2);
            }
        }
    }

    std::vector<Node> nodes;
    int root = -1;
    uint64_t relevant = 0;      // atoms that occur in some Inf/Fin
    long budget = long(1) << 24;
};

class ZielonkaTree final {
public:
    struct Node {
        uint64_t atoms;
        bool accepting;
        int parent;
        int depth;
        std::vector<int> children;
    };

    ZielonkaTree(AtomCondition& condition, uint64_t all)
    {
        nodes.push_back(Node{all, condition.accepts(all), -1, 0, {}});
        for (size_t n = 0; n < nodes.size(); n++) {
            const uint64_t atoms = nodes[n].atoms;
            const bool accepting = nodes[n].accepting;
            // the maximal subsets with the other acceptance
            std::vector<uint64_t> children;
            uint64_t x;
            while (condition.find(atoms, 0, !accepting, children, x)) {
                std::vector<uint64_t> larger{x};
                uint64_t y;
                while (condition.find(atoms & ~x, x, !accepting, larger, y)) larger[0] = x = y;
                children.push_back(x);
            }
            for (uint64_t c : children) {
                if (nodes.size() == MAX_NODES) {
                    conversionError("the Zielonka tree of the acceptance condition is too large");
                }
                nodes[n].children.push_back(static_cast<int>(nodes.size()));
                nodes.push_back(Node{c, !accepting, static_cast<int>(n), nodes[n].depth + 1, {}});
            }
        }
        int depth = 0;
        for (const Node& n : nodes) {
            depth = std::max(depth, n.depth);
            if (n.children.empty()) {
                leafOf.push_back(static_cast<int>(&n - nodes.data()));
            }
        }
        // Node n gets priority maxDepth - depth(n) + base: the shallower, the
        // higher, and base makes the priorities of accepting nodes even.
        base = (depth % 2 == 0) == nodes[0].accepting ? 0 : 1;
        maxDepth = depth;
        leafIndex.assign(nodes.size(), -1);
        for (size_t i = 0; i < leafOf.size(); i++) leafIndex[leafOf[i]] = static_cast<int>(i);
    }

    size_t leaves() const { return leafOf.size(); }

    // The initial state (leaf) of the Zielonka tree automaton.
    int initial() const { return leftmostLeaf(0); }

    // The priority and successor leaf when the Zielonka tree automaton reads
    // the atom in the given leaf: for the deepest ancestor n of the leaf that
    // contains the atom, the priority of n and the leftmost leaf below the
    // next child of n (cyclically) after the one towards the leaf.
    std::pair<int, int> step(int leaf, int atom) const
    {
        int n = leafOf[leaf];
        int child = -1;
        while ((nodes[n].atoms >> atom & 1) == 0) {
            child = n;
            n = nodes[n].parent;
        }
        const int priority = maxDepth - nodes[n].depth + base;
        if (child == -1) return {priority, leaf};
        const std::vector<int>& cs = nodes[n].children;
        const size_t i = static_cast<size_t>(std::find(cs.begin(), cs.end(), child) - cs.begin());
        return {priority, leftmostLeaf(cs[(i + 1) % cs.size()])};
    }

private:
    static constexpr size_t MAX_NODES = size_t(1) << 14;

    int leftmostLeaf(int n) const
    {
        while (!nodes[n].children.empty()) n = nodes[n].children[0];
        return leafIndex[n];
    }

    std::vector<Node> nodes;
    std::vector<int> leafOf;        // leaf number -> node
    std::vector<int> leafIndex;     // node -> leaf number (or -1)
    int maxDepth = 0;
    int base = 0;
};

uint64_t setMask(const int* sets, int count)
{
    uint64_t m = 0;
    for (int i = 0; i < count; i++) m |= uint64_t(1) << sets[i];
    return m;
}

// The canonical "parity max even" condition over priorities 0..p (malloc'd).
TreePtr parityMaxEven(int p)
{
    const bool good = p % 2 == 0;
    TreePtr term = mallocNode(good ? NT_INF : NT_FIN, mallocNode(NT_SET, nullptr, nullptr, p), nullptr, -1);
    if (p == 0) return term;
    TreePtr rest = parityMaxEven(p - 1);
    return mallocNode(good ? NT_OR : NT_AND, std::move(term), std::move(rest), -1);
}

// Renumber priorities in place, keeping their order and parities, so that
// none are wasted; returns the number of priorities.
int compactPriorities(std::vector<int>& priorities)
{
    std::vector<int> used(priorities);
    std::sort(used.begin(), used.end());
    used.erase(std::unique(used.begin(), used.end()), used.end());
    if (used.empty()) return 1;
    std::vector<int> renumbered(used.back() + 1, 0);
    int n = used[0] % 2;
    renumbered[used[0]] = n;
    for (size_t i = 1; i < used.size(); i++) {
        if (used[i] % 2 != used[i - 1] % 2) n++;
        renumbered[used[i]] = n;
    }
    for (int& p : priorities) p = renumbered[p];
    return n + 1;
}

// The edges of the automaton in order, with their atoms.
struct Atoms {
    std::vector<uint64_t> masks;    // atom -> its acceptance sets
    std::vector<int> ofEdge;        // edge (in state order) -> atom
};

Atoms collectAtoms(const HoaData* data)
{
    Atoms a;
    std::unordered_map<uint64_t, int> index;
    for (int i = 0; i < data->noStates; i++) {
        const State& s = data->states[i];
        const uint64_t stateSets = setMask(s.accSig, s.noAccSig);
        for (int j = 0; j < s.noTrans; j++) {
            const Transition& t = s.transitions[j];
            const uint64_t sets = stateSets | setMask(t.accSig, t.noAccSig);
            const auto inserted = index.emplace(sets, static_cast<int>(a.masks.size()));
            if (inserted.second) a.masks.push_back(sets);
            a.ofEdge.push_back(inserted.first->second);
        }
    }
    if (a.masks.size() > 64) conversionError("more than 64 distinct combinations of acceptance sets on edges");
    return a;
}

// Set the acceptance condition to "parity max even" over noSets priorities
// and update the properties: acceptance is on edges, one set per edge.
void setParityAcceptance(HoaData* data, int noSets, bool product)
{
    std::vector<std::string> props;
    for (const std::string& p : properties(data)) {
        // a product may break very-weakness and tightness
        if (p != "state-acc" && !(product && (p == "very-weak" || p == "tight"))) props.push_back(p);
    }
    addProperty(props, "trans-acc");
    addProperty(props, "colored");

    TreePtr acc = parityMaxEven(noSets - 1);
    char* name = mallocString("parity");
    char** parameters = nullptr;
    try {
        parameters = mallocStrings({"max", "even", std::to_string(noSets)});
        replaceStrings(data->properties, data->noProps, props);
    } catch (...) {
        if (parameters != nullptr) {
            for (int i = 0; i < 3; i++) free(parameters[i]);
            free(parameters);
        }
        free(name);
        throw;
    }
    freeTree(data->acc);
    data->acc = acc.release();
    data->noAccSets = noSets;
    free(data->accNameID);
    data->accNameID = name;
    for (int i = 0; i < data->noANPs; i++) free(data->accNameParameters[i]);
    free(data->accNameParameters);
    data->accNameParameters = parameters;
    data->noANPs = 3;
}

// Recolor the edges in place: edge e gets priority[atom of e].
void recolor(HoaData* data, const Atoms& atoms, std::vector<int> priority)
{
    const int noSets = compactPriorities(priority);
    MallocList fresh;
    std::vector<int*> sigs;
    for (int i = 0; i < data->noStates; i++) {
        const State& s = data->states[i];
        for (int j = 0; j < s.noTrans; j++) {
            if (s.transitions[j].noAccSig == 0) sigs.push_back(fresh.allocate<int>(1));
        }
    }
    setParityAcceptance(data, noSets, false);
    fresh.release();
    size_t e = 0;
    size_t next = 0;
    for (int i = 0; i < data->noStates; i++) {
        State& s = data->states[i];
        for (int j = 0; j < s.noTrans; j++, e++) {
            Transition& t = s.transitions[j];
            if (t.noAccSig == 0) t.accSig = sigs[next++];
            t.accSig[0] = priority[atoms.ofEdge[e]];
            t.noAccSig = 1;
        }
        free(s.accSig);
        s.accSig = nullptr;
        s.noAccSig = 0;
    }
}

// Replace the automaton by its product with the Zielonka tree automaton
// (only the reachable part).
void multiply(HoaData* data, const Atoms& atoms, const ZielonkaTree& tree)
{
    Arena* arena = static_cast<Arena*>(data->arena);
    if (arena == nullptr) conversionError("adding memory requires an automaton read by the parser");
    const size_t leaves = tree.leaves();
    if (leaves * atoms.masks.size() > (size_t(1) << 22)) conversionError("the Zielonka tree has too many leaves");

    std::vector<size_t> firstEdge(data->noStates + 1, 0);
    for (int i = 0; i < data->noStates; i++) {
        const State& s = data->states[i];
        for (int j = 0; j < s.noTrans; j++) {
            if (s.transitions[j].noSucc != 1) conversionError("adding memory to alternating automata is not supported");
        }
        firstEdge[i + 1] = firstEdge[i] + static_cast<size_t>(s.noTrans);
    }

    // breadth-first over the reachable pairs (state, leaf)
    std::vector<std::pair<int, int>> pairs;
    std::unordered_map<uint64_t, int> id;
    auto pairId = [&](int state, int leaf) {
        const uint64_t key = static_cast<uint64_t>(state) * leaves + static_cast<uint64_t>(leaf);
        if (pairs.size() == INT_MAX) conversionError("the product is too large");
        const auto inserted = id.emplace(key, static_cast<int>(pairs.size()));
        if (inserted.second) pairs.emplace_back(state, leaf);
        return inserted.first->second;
    };
    std::vector<int> start;
    for (int i = 0; i < data->noStart; i++) start.push_back(pairId(data->start[i], tree.initial()));
    std::vector<std::pair<int, int>> steps(leaves * atoms.masks.size(), {-1, -1});
    std::vector<int> successor;
    std::vector<int> priority;
    for (size_t i = 0; i < pairs.size(); i++) {
        const int q = pairs[i].first;
        const int leaf = pairs[i].second;
        const State& s = data->states[q];
        for (int j = 0; j < s.noTrans; j++) {
            const int atom = atoms.ofEdge[firstEdge[q] + static_cast<size_t>(j)];
            std::pair<int, int>& st = steps[static_cast<size_t>(leaf) * atoms.masks.size() + atom];
            if (st.first < 0) st = tree.step(leaf, atom);
            priority.push_back(st.first);
            successor.push_back(pairId(s.transitions[j].successors[0], st.second));
        }
    }
    const int noSets = compactPriorities(priority);

    // Build the new states; edges share the labels of the original edges.
    MallocList fresh;
    State* states = fresh.allocate<State>(pairs.size());
    int* initial = fresh.allocate<int>(start.size());
    std::copy(start.begin(), start.end(), initial);
    size_t e = 0;
    for (size_t i = 0; i < pairs.size(); i++) {
        const State& original = data->states[pairs[i].first];
        State& s = states[i];
        s.id = static_cast<int>(i);
        s.name = nullptr;
        s.label = original.label;
        s.accSig = nullptr;
        s.noAccSig = 0;
        s.noTrans = original.noTrans;
        s.transitions = arena->allocate<Transition>(static_cast<size_t>(s.noTrans));
        for (int j = 0; j < s.noTrans; j++, e++) {
            Transition& t = s.transitions[j];
            t.label = original.transitions[j].label;
            t.successors = arena->allocate<int>(1);
            t.successors[0] = successor[e];
            t.noSucc = 1;
            t.accSig = fresh.allocate<int>(1);
            t.accSig[0] = priority[e];
            t.noAccSig = 1;
        }
    }
    setParityAcceptance(data, noSets, true);
    fresh.release();

    for (int i = 0; i < data->noStates; i++) {
        State& s = data->states[i];
        free(s.accSig);
        for (int j = 0; j < s.noTrans; j++) free(s.transitions[j].accSig);
    }
    free(data->states);
    free(data->start);
    data->states = states;
    data->noStates = static_cast<int>(pairs.size());
    data->start = start.empty() ? nullptr : initial;
    if (start.empty()) free(initial);
    data->noStart = static_cast<int>(start.size());
}

void paritize(HoaData* data, bool addMemory)
{
    bool max;
    short good;
    if (parityName(data, &max, &good) == 0
        && isParityCondition(data->acc, data->noAccSets, max, good, max ? data->noAccSets - 1 : 0, true)
        && coloring(data) != Coloring::NONE) {
        // already usable as it is
        std::vector<std::string> props = properties(data);
        if (std::find(props.begin(), props.end(), "colored") == props.end()) {
            addProperty(props, "colored");
            replaceStrings(data->properties, data->noProps, props);
        }
        return;
    }
    if (data->acc == nullptr) conversionError("there is no acceptance condition");
    if (data->noAccSets > 64) conversionError("more than 64 acceptance sets");

    const Atoms atoms = collectAtoms(data);
    AtomCondition condition;
    if (!condition.build(data->acc, atoms.masks)) conversionError("unsupported acceptance condition");
    const uint64_t all = atoms.masks.size() == 64 ? ~uint64_t(0) : (uint64_t(1) << atoms.masks.size()) - 1;
    const ZielonkaTree tree(condition, all);

    if (tree.leaves() == 1) {
        // a chain: priorities depend only on the atom
        std::vector<int> priority;
        for (size_t a = 0; a < atoms.masks.size(); a++) priority.push_back(tree.step(0, static_cast<int>(a)).first);
        recolor(data, atoms, std::move(priority));
    } else if (addMemory) {
        multiply(data, atoms, tree);
    } else {
        conversionError("the acceptance condition is not a parity condition on this automaton; it needs memory "
                        "(its Zielonka tree has " + std::to_string(tree.leaves()) + " leaves)", 104);
    }
}

}  // namespace

extern "C" int paritizeHoa(HoaData* data, bool addMemory)
{
    try {
        paritize(data, addMemory);
        return 0;
    } catch (const ConversionError& e) {
        fprintf(stderr, "Unsupported acceptance condition: %s\n", e.message.c_str());
        return e.code;
    } catch (const std::bad_alloc&) {
        fprintf(stderr, "Out of memory\n");
        return 100;
    }
}
