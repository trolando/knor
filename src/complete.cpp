/**
 * Copyright Tom van Dijk
 */

/*
 * Completion of automata that do not declare "deterministic" and "complete"
 * (completeHoa, see complete.hpp). The labels of a state are compared as
 * BDDs: its edges with the same successor and the same acceptance sets form
 * a group, the automaton is deterministic iff the groups of every state are
 * disjoint, and complete iff they cover every valuation.
 *
 * A valuation without an edge rejects the run, so an incomplete automaton
 * gets a rejecting sink state, which the incomplete states reach with the
 * missing valuations:
 *  - a parity automaton with one acceptance set per edge (or per state) keeps
 *    a parity condition of the same kind: the sink gets the most significant
 *    priority if that one is losing, and otherwise a new losing priority
 *    above all others (max parity), or the priorities move up by one and the
 *    sink gets the new priority 0 (min parity);
 *  - any other automaton gets a fresh acceptance set s on the sink's loop and
 *    the condition "acc & Fin(s)", which paritizeHoa turns into parity.
 */

#include <algorithm>
#include <cstdio>
#include <deque>
#include <map>
#include <new>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <sylvan.h>

#include "complete.hpp"
#include "hoainternal.hpp"

using namespace hoa;
using namespace sylvan;

namespace {

struct CompletionError {
    int code;
    std::string message;
};

[[noreturn]] void completionError(std::string message, int code = 100)
{
    throw CompletionError{code, std::move(message)};
}

// Protected BDD slots with stable addresses.
class Slots final {
public:
    Slots() = default;
    Slots(const Slots&) = delete;
    Slots& operator=(const Slots&) = delete;
    ~Slots()
    {
        for (MTBDD& s : slots) mtbdd_unprotect(&s);
    }

    MTBDD& operator[](size_t i)
    {
        while (slots.size() <= i) {
            slots.push_back(mtbdd_false);
            mtbdd_protect(&slots.back());
        }
        return slots[i];
    }

private:
    std::deque<MTBDD> slots;
};

// A local on Sylvan's reference stack for the lifetime of the guard.
class Ref final {
public:
    explicit Ref(MTBDD& r) { mtbdd_refs_pushptr(&r); }
    ~Ref() { mtbdd_refs_popptr(1); }
    Ref(const Ref&) = delete;
    Ref& operator=(const Ref&) = delete;
};

// BDDs of label expressions, with proposition i as variable i.
class Labels final {
public:
    explicit Labels(const HoaData* data) : noAPs(data->noAPs)
    {
        for (int i = 0; i < data->noAliases; i++) aliases.emplace(data->aliases[i].alias, data->aliases[i].labelExpr);
    }

    // The BDD of t (true if t is NULL) into result, which the caller protects.
    void bdd(const BTree* t, MTBDD& result)
    {
        if (t == nullptr) {
            result = mtbdd_true;
            return;
        }
        switch (t->type) {
            case NT_BOOL:
                result = t->id != 0 ? mtbdd_true : mtbdd_false;
                return;
            case NT_AP:
                if (t->id < 0 || t->id >= noAPs) completionError("label with an undeclared proposition");
                result = sylvan_ithvar(static_cast<BDDVAR>(t->id));
                return;
            case NT_NOT:
                bdd(t->left, result);
                result = sylvan_not(result);
                return;
            case NT_ALIAS:
                result = alias(t->alias);
                return;
            case NT_AND:
            case NT_OR:
                chain(t, result);
                return;
            default:
                completionError("unexpected operator in a label");
        }
    }

private:
    // A chain of the same operator, evaluated without recursion along it.
    void chain(const BTree* t, MTBDD& result)
    {
        const bool conjunction = t->type == NT_AND;
        std::vector<const BTree*> operands;
        std::vector<const BTree*> open{t};
        while (!open.empty()) {
            const BTree* n = open.back();
            open.pop_back();
            if (n->type == t->type) {
                open.push_back(n->right);
                open.push_back(n->left);
            } else {
                operands.push_back(n);
            }
        }
        const MTBDD dominant = conjunction ? mtbdd_false : mtbdd_true;
        MTBDD value = conjunction ? mtbdd_true : mtbdd_false;
        MTBDD operand = mtbdd_false;
        Ref value_ref(value);
        Ref operand_ref(operand);
        for (const BTree* o : operands) {
            bdd(o, operand);
            value = conjunction ? sylvan_and(value, operand) : sylvan_or(value, operand);
            if (value == dominant) break;
        }
        result = value;
    }

    MTBDD alias(const char* name)
    {
        const auto known = aliasSlot.find(name);
        if (known != aliasSlot.end()) return aliasBdds[known->second];
        const auto a = aliases.find(name);
        if (a == aliases.end()) completionError(std::string("undefined alias @") + name);
        MTBDD value = mtbdd_false;
        Ref value_ref(value);
        bdd(a->second, value);
        const size_t i = aliasSlot.size();
        aliasBdds[i] = value;
        aliasSlot.emplace(name, i);
        return value;
    }

    int noAPs;
    std::unordered_map<std::string, const BTree*> aliases;
    std::unordered_map<std::string, size_t> aliasSlot;
    Slots aliasBdds;
};

// The canonical parity condition with n priorities (as isParityCondition, strict).
TreePtr canonicalParity(bool max, short good, int n)
{
    auto term = [good](int p) {
        return mallocNode(p % 2 == good ? NT_INF : NT_FIN, mallocNode(NT_SET, nullptr, nullptr, p), nullptr, -1);
    };
    TreePtr condition = term(max ? 0 : n - 1);
    for (int k = 1; k < n; k++) {
        const int p = max ? k : n - 1 - k;
        condition = mallocNode(p % 2 == good ? NT_OR : NT_AND, term(p), std::move(condition), -1);
    }
    return condition;
}

BTree* arenaNode(Arena* arena, NodeType type, BTree* left, BTree* right, int id)
{
    BTree* n = arena->allocate<BTree>(1);
    n->left = left;
    n->right = right;
    n->alias = nullptr;
    n->id = id;
    n->type = type;
    return n;
}

// Add the rejecting sink (see the top of this file).
void addSink(HoaData* data, Arena* arena, const std::vector<int>& incomplete)
{
    if (data->acc == nullptr) completionError("there is no acceptance condition");
    bool max = true;
    short good = 0;
    const Coloring colors = coloring(data);
    const bool parity = data->noAccSets >= 1 && colors != Coloring::NONE && parityName(data, &max, &good) == 0
                        && isParityCondition(data->acc, data->noAccSets, max, good, max ? data->noAccSets - 1 : 0, true);
    int noSets = data->noAccSets;
    int sinkSet = noSets;
    bool shift = false;
    short newGood = good;
    if (!parity) {
        noSets++;
    } else if ((max ? noSets - 1 : 0) % 2 != good) {
        sinkSet = max ? noSets - 1 : 0;
    } else if (max) {
        noSets++;
    } else {
        // all priorities move up by one, which flips their parity
        shift = true;
        newGood = static_cast<short>(1 - good);
        sinkSet = 0;
        noSets++;
    }
    const bool onEdges = !parity || colors == Coloring::EDGES;

    MallocList mallocd;
    const int sink = data->noStates;
    State* states = mallocd.allocate<State>(static_cast<size_t>(data->noStates) + 1);
    std::copy(data->states, data->states + data->noStates, states);
    std::vector<int*> sinkSigs;
    for (const int i : incomplete) {
        State& s = states[i];
        Transition* transitions = arena->allocate<Transition>(static_cast<size_t>(s.noTrans) + 1);
        std::copy(s.transitions, s.transitions + s.noTrans, transitions);
        BTree* covered = nullptr;
        for (int j = 0; j < s.noTrans; j++) {
            if (transitions[j].label == nullptr) transitions[j].label = s.label;
            covered = j == 0 ? transitions[j].label : arenaNode(arena, NT_OR, covered, transitions[j].label, -1);
        }
        Transition& missing = transitions[s.noTrans];
        missing.label = covered == nullptr ? arenaNode(arena, NT_BOOL, nullptr, nullptr, 1)
                                           : arenaNode(arena, NT_NOT, covered, nullptr, -1);
        missing.successors = arena->allocate<int>(1);
        missing.successors[0] = sink;
        missing.noSucc = 1;
        // a parity automaton colored on its edges needs a priority on every edge
        missing.noAccSig = parity && colors == Coloring::EDGES ? 1 : 0;
        missing.accSig = missing.noAccSig != 0 ? mallocd.allocate<int>(1) : nullptr;
        if (missing.accSig != nullptr) sinkSigs.push_back(missing.accSig);
        s.label = nullptr;
        s.transitions = transitions;
        s.noTrans++;
    }
    State& s = states[sink];
    s.id = sink;
    s.name = nullptr;
    s.label = nullptr;
    s.noAccSig = onEdges ? 0 : 1;
    s.accSig = onEdges ? nullptr : mallocd.allocate<int>(1);
    if (s.accSig != nullptr) sinkSigs.push_back(s.accSig);
    s.noTrans = 1;
    s.transitions = arena->allocate<Transition>(1);
    Transition& loop = s.transitions[0];
    loop.label = arenaNode(arena, NT_BOOL, nullptr, nullptr, 1);
    loop.successors = arena->allocate<int>(1);
    loop.successors[0] = sink;
    loop.noSucc = 1;
    loop.noAccSig = onEdges ? 1 : 0;
    loop.accSig = onEdges ? mallocd.allocate<int>(1) : nullptr;
    if (loop.accSig != nullptr) sinkSigs.push_back(loop.accSig);

    TreePtr acc;
    char** parameters = nullptr;
    if (!parity) {
        acc = mallocNode(NT_AND, nullptr, mallocNode(NT_FIN, mallocNode(NT_SET, nullptr, nullptr, sinkSet), nullptr, -1), -1);
    } else if (noSets != data->noAccSets) {
        acc = canonicalParity(max, newGood, noSets);
        parameters = mallocStrings({max ? "max" : "min", newGood != 0 ? "odd" : "even", std::to_string(noSets)});
    }

    // the rest cannot fail
    mallocd.release();
    if (shift) {
        for (int i = 0; i < data->noStates; i++) {
            const State& old = data->states[i];
            for (int k = 0; k < old.noAccSig; k++) old.accSig[k]++;
            for (int j = 0; j < old.noTrans; j++) {
                for (int k = 0; k < old.transitions[j].noAccSig; k++) old.transitions[j].accSig[k]++;
            }
        }
    }
    for (int* sig : sinkSigs) sig[0] = sinkSet;
    if (!parity) {
        acc->left = data->acc;
        data->acc = acc.release();
        free(data->accNameID);
        data->accNameID = nullptr;
        for (int i = 0; i < data->noANPs; i++) free(data->accNameParameters[i]);
        free(data->accNameParameters);
        data->accNameParameters = nullptr;
        data->noANPs = 0;
    } else if (acc != nullptr) {
        freeTree(data->acc);
        data->acc = acc.release();
        for (int i = 0; i < data->noANPs; i++) free(data->accNameParameters[i]);
        free(data->accNameParameters);
        data->accNameParameters = parameters;
        data->noANPs = 3;
    }
    data->noAccSets = noSets;
    free(data->states);
    data->states = states;
    data->noStates++;
}

void complete(HoaData* data)
{
    std::vector<std::string> props = properties(data);
    const bool deterministic = std::find(props.begin(), props.end(), "deterministic") != props.end();
    const bool complete = std::find(props.begin(), props.end(), "complete") != props.end();
    if (deterministic && complete) return;
    Arena* arena = static_cast<Arena*>(data->arena);
    if (arena == nullptr) completionError("checking the automaton requires an automaton read by the parser");

    Labels labels(data);
    Slots slots;  // 0: covered, 1: scratch, 2...: the groups of the state
    std::vector<int> incomplete;
    std::map<std::pair<int, std::vector<int>>, size_t> group;
    for (int i = 0; i < data->noStates; i++) {
        const State& s = data->states[i];
        group.clear();
        for (int j = 0; j < s.noTrans; j++) {
            const Transition& t = s.transitions[j];
            if (t.noSucc != 1) completionError("alternating automata are not supported");
            std::vector<int> sets(t.accSig, t.accSig + t.noAccSig);
            std::sort(sets.begin(), sets.end());
            sets.erase(std::unique(sets.begin(), sets.end()), sets.end());
            const auto g = group.emplace(std::make_pair(t.successors[0], std::move(sets)), group.size());
            MTBDD& edges = slots[2 + g.first->second];
            if (g.second) edges = mtbdd_false;
            labels.bdd(t.label != nullptr ? t.label : s.label, slots[1]);
            edges = sylvan_or(edges, slots[1]);
        }
        MTBDD& covered = slots[0];
        covered = mtbdd_false;
        for (size_t g = 0; g < group.size(); g++) {
            if (!deterministic) {
                slots[1] = sylvan_and(slots[2 + g], covered);
                if (slots[1] != mtbdd_false) {
                    completionError("the automaton is not deterministic: state " + std::to_string(i)
                                    + " has two edges for the same valuation", 200);
                }
            }
            covered = sylvan_or(covered, slots[2 + g]);
        }
        if (!complete && covered != mtbdd_true) incomplete.push_back(i);
    }
    addProperty(props, "deterministic");
    addProperty(props, "complete");
    if (!incomplete.empty()) addSink(data, arena, incomplete);
    replaceStrings(data->properties, data->noProps, props);
}

}  // namespace

int completeHoa(HoaData* data)
{
    try {
        complete(data);
        return 0;
    } catch (const CompletionError& e) {
        fprintf(stderr, "Unsupported automaton: %s\n", e.message.c_str());
        return e.code;
    } catch (const std::bad_alloc&) {
        fprintf(stderr, "Out of memory\n");
        return 100;
    }
}
