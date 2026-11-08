/**
 * Copyright Tom van Dijk
 */

/*
 * aigcheck: check an AIGER controller against the eHOA parity automaton it was
 * synthesized from.
 *
 *     aigcheck automaton.ehoa controller.aag
 *
 * The controller is a SYNTCOMP Mealy machine: one AIGER input per
 * uncontrollable AP and one AIGER output per controllable AP, matched by
 * symbol name, latches initialised to 0 or 1, outputs computed from the
 * current inputs and latches. The automaton must be what knor accepts
 * (isParityGFG): a complete deterministic colored parity automaton with a
 * unique initial state, with the priority on the edge or on the target state.
 *
 * The checker explores the product of the controller's latch valuations with
 * the automaton states over all input valuations and reports FAIL iff some
 * reachable cycle has a dominant priority (max or min according to acc-name)
 * of the wrong parity, i.e. iff the controller has a trace that the automaton
 * rejects. Input valuations are enumerated as cubes: inputs are assigned one
 * at a time while the circuit and the labels are evaluated in three-valued
 * logic, and a cube is closed as soon as the edge taken and the next latch
 * values are determined. This is exact and independent of knor's BDD and AIG
 * code. The properties "deterministic" and "complete" are checked on every
 * reached valuation: one without an edge is reported as FAIL (the automaton is
 * not complete, and knor's game has no move there), one with several edges is
 * an error.
 *
 * Exit codes: 0 PASS, 1 FAIL, 2 error (usage, parsing, malformed input).
 */

#include <algorithm>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <pthread.h>

extern "C" {
#include "aiger.h"
}
#include "simplehoa.h"

namespace {

[[noreturn]] void fail(const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    fputs("aigcheck: ", stderr);
    vfprintf(stderr, fmt, args);
    fputc('\n', stderr);
    va_end(args);
    exit(2);
}

// Three-valued logic: 0, 1 and X (unknown).
const uint8_t VX = 2;

uint8_t not3(uint8_t a) { return a == VX ? VX : static_cast<uint8_t>(1 - a); }
uint8_t and3(uint8_t a, uint8_t b) { return (a == 0 || b == 0) ? 0 : (a == 1 && b == 1) ? 1 : VX; }
uint8_t or3(uint8_t a, uint8_t b) { return not3(and3(not3(a), not3(b))); }

/**
 * The specification: the parsed automaton, its parity condition, and for
 * every edge its priority and the APs in its label.
 */
class Spec final {
public:
    explicit Spec(HoaData& data) : data_(data)
    {
        // As knor does: rewrite other acceptance conditions (Buchi, Rabin, ...)
        // into canonical parity; this preserves the language.
        if (paritizeHoa(&data, true) != 0) exit(2);
        if (isParityGFG(&data, &is_max_, &good_) != 0) exit(2);
        if (data.noAccSets < 1) fail("the automaton has no acceptance sets");
        for (int i = 0; i < data.noAliases; i++) aliases_[data.aliases[i].alias] = data.aliases[i].labelExpr;
        edge_prio_.resize(data.noStates);
        edge_aps_.resize(data.noStates);
        for (int s = 0; s < data.noStates; s++) {
            const State& st = data.states[s];
            for (int j = 0; j < st.noTrans; j++) {
                const Transition& t = st.transitions[j];
                if (t.noSucc != 1) fail("state %d has an edge without exactly one successor", s);
                if (label(s, j) == nullptr) fail("state %d has an edge without a label", s);
                int p;
                if (t.noAccSig == 1) {
                    p = t.accSig[0];
                } else if (t.noAccSig == 0 && data.states[t.successors[0]].noAccSig == 1) {
                    p = data.states[t.successors[0]].accSig[0];
                } else {
                    fail("edge %d of state %d is not colored by exactly one acceptance set (on the edge or on the target state)", j, s);
                }
                if (p < 0 || p >= data.noAccSets) fail("acceptance set %d of state %d is out of range", p, s);
                edge_prio_[s].push_back(p);
                std::vector<bool> used(data.noAPs, false);
                collectAPs(label(s, j), used);
                std::vector<int> aps;
                for (int a = 0; a < data.noAPs; a++) if (used[a]) aps.push_back(a);
                edge_aps_[s].push_back(std::move(aps));
            }
        }
    }

    const HoaData& data() const { return data_; }
    bool isMax() const { return is_max_; }
    bool isGood(int p) const { return p % 2 == good_; }
    int start() const { return data_.start[0]; }
    int edges(int s) const { return data_.states[s].noTrans; }
    int target(int s, int j) const { return data_.states[s].transitions[j].successors[0]; }
    int prio(int s, int j) const { return edge_prio_[s][j]; }
    const std::vector<int>& aps(int s, int j) const { return edge_aps_[s][j]; }

    const BTree* label(int s, int j) const
    {
        const State& st = data_.states[s];
        return st.label != nullptr ? st.label : st.transitions[j].label;
    }

    // A (partial) valuation of the APs as "a !b ?c" (? for unassigned).
    std::string valuationString(const std::vector<uint8_t>& ap) const
    {
        std::string s;
        for (int a = 0; a < data_.noAPs; a++) {
            if (!s.empty()) s += ' ';
            s += ap[a] == 1 ? "" : ap[a] == 0 ? "!" : "?";
            s += data_.aps[a];
        }
        return s.empty() ? "(no APs)" : s;
    }

    uint8_t eval(const BTree* t, const std::vector<uint8_t>& ap) const
    {
        switch (t->type) {
            case NT_BOOL: return t->id != 0 ? 1 : 0;
            case NT_AP: return ap[checkedAP(t->id)];
            case NT_NOT: return not3(eval(t->left, ap));
            case NT_AND: return and3(eval(t->left, ap), eval(t->right, ap));
            case NT_OR: return or3(eval(t->left, ap), eval(t->right, ap));
            case NT_ALIAS: return eval(alias(t->alias), ap);
            default: fail("unsupported operator in a label");
        }
    }

private:
    int checkedAP(int id) const
    {
        if (id < 0 || id >= data_.noAPs) fail("atomic proposition %d in a label is out of range", id);
        return id;
    }

    const BTree* alias(const char* name) const
    {
        auto it = aliases_.find(name);
        if (it == aliases_.end()) fail("unknown alias @%s", name);
        return it->second;
    }

    void collectAPs(const BTree* t, std::vector<bool>& used) const
    {
        switch (t->type) {
            case NT_BOOL: return;
            case NT_AP: used[checkedAP(t->id)] = true; return;
            case NT_NOT: collectAPs(t->left, used); return;
            case NT_AND:
            case NT_OR: collectAPs(t->left, used); collectAPs(t->right, used); return;
            case NT_ALIAS: collectAPs(alias(t->alias), used); return;
            default: fail("unsupported operator in a label");
        }
    }

    const HoaData& data_;
    bool is_max_ = true;
    short good_ = 0;
    std::map<std::string, const BTree*> aliases_;
    std::vector<std::vector<int>> edge_prio_;
    std::vector<std::vector<std::vector<int>>> edge_aps_;
};

/**
 * The controller: the AIG, the bijections between its inputs/outputs and the
 * uncontrollable/controllable APs, and per AIG variable the set of inputs in
 * its structural cone (as a bitset over input indices).
 */
class Controller final {
public:
    Controller(const char* filename, const HoaData& data) : aig_(aiger_init())
    {
        if (const char* err = aiger_open_and_read_from_file(aig_, filename)) fail("%s", err);
        if (aig_->num_bad + aig_->num_constraints + aig_->num_justice + aig_->num_fairness > 0) {
            fail("the controller has bad/constraint/justice/fairness sections");
        }
        aiger_reencode(aig_);   // inputs, latches, then AND gates in topological order
        for (unsigned i = 0; i < aig_->num_latches; i++) {
            const unsigned reset = aig_->latches[i].reset;
            if (reset != 0 && reset != 1) fail("latch %u is not initialised to 0 or 1", aig_->latches[i].lit);
        }

        std::map<std::string, int> ap_index;
        for (int i = 0; i < data.noAPs; i++) ap_index[data.aps[i]] = i;
        std::vector<bool> controllable(data.noAPs, false);
        for (int i = 0; i < data.noCntAPs; i++) controllable[data.cntAPs[i]] = true;
        ap_input_.assign(data.noAPs, -1);
        ap_output_.assign(data.noAPs, -1);
        auto bind = [&](const aiger_symbol& sym, unsigned idx, bool is_output) {
            const char* kind = is_output ? "output" : "input";
            if (sym.name == nullptr) fail("%s %u has no name", kind, idx);
            auto it = ap_index.find(sym.name);
            if (it == ap_index.end()) fail("%s \"%s\" is not an atomic proposition", kind, sym.name);
            const int ap = it->second;
            if (controllable[ap] != is_output) fail("%s \"%s\" is %s", kind, sym.name, is_output ? "not a controllable AP" : "a controllable AP");
            std::vector<int>& slot = is_output ? ap_output_ : ap_input_;
            if (slot[ap] != -1) fail("%s \"%s\" occurs twice", kind, sym.name);
            slot[ap] = static_cast<int>(idx);
        };
        for (unsigned i = 0; i < aig_->num_inputs; i++) bind(aig_->inputs[i], i, false);
        for (unsigned i = 0; i < aig_->num_outputs; i++) bind(aig_->outputs[i], i, true);
        for (int a = 0; a < data.noAPs; a++) {
            if (!controllable[a] && ap_input_[a] == -1) fail("uncontrollable AP \"%s\" has no input", data.aps[a]);
            if (controllable[a] && ap_output_[a] == -1) fail("controllable AP \"%s\" has no output", data.aps[a]);
        }
        in_ap_.resize(aig_->num_inputs);
        out_ap_.resize(aig_->num_outputs);
        for (int a = 0; a < data.noAPs; a++) {
            if (ap_input_[a] != -1) in_ap_[ap_input_[a]] = a;
            if (ap_output_[a] != -1) out_ap_[ap_output_[a]] = a;
        }

        words_ = (aig_->num_inputs + 63) / 64;
        support_.assign(aig_->maxvar + 1, std::vector<uint64_t>(words_, 0));
        for (unsigned i = 0; i < aig_->num_inputs; i++) support_[aig_->inputs[i].lit >> 1][i / 64] |= uint64_t(1) << (i % 64);
        for (unsigned i = 0; i < aig_->num_ands; i++) {
            const aiger_and& g = aig_->ands[i];
            for (unsigned w = 0; w < words_; w++) {
                support_[g.lhs >> 1][w] = support_[g.rhs0 >> 1][w] | support_[g.rhs1 >> 1][w];
            }
        }
    }

    ~Controller() { aiger_reset(aig_); }
    Controller(const Controller&) = delete;
    Controller& operator=(const Controller&) = delete;

    unsigned inputs() const { return aig_->num_inputs; }
    unsigned outputs() const { return aig_->num_outputs; }
    unsigned latches() const { return aig_->num_latches; }
    unsigned words() const { return words_; }
    int inputAP(unsigned i) const { return in_ap_[i]; }
    int outputAP(unsigned i) const { return out_ap_[i]; }
    int apInput(int ap) const { return ap_input_[ap]; }
    int apOutput(int ap) const { return ap_output_[ap]; }
    std::vector<uint8_t> initialLatches() const
    {
        std::vector<uint8_t> init(latches());
        for (unsigned i = 0; i < latches(); i++) init[i] = static_cast<uint8_t>(aig_->latches[i].reset);
        return init;
    }

    // Three-valued simulation of one step: the value of every AIG variable.
    std::vector<uint8_t> simulate(const std::vector<uint8_t>& input, const std::vector<uint8_t>& latch) const
    {
        std::vector<uint8_t> v(aig_->maxvar + 1, 0);
        for (unsigned i = 0; i < aig_->num_inputs; i++) v[aig_->inputs[i].lit >> 1] = input[i];
        for (unsigned i = 0; i < aig_->num_latches; i++) v[aig_->latches[i].lit >> 1] = latch[i];
        for (unsigned i = 0; i < aig_->num_ands; i++) {
            const aiger_and& g = aig_->ands[i];
            v[g.lhs >> 1] = and3(value(v, g.rhs0), value(v, g.rhs1));
        }
        return v;
    }

    uint8_t output(const std::vector<uint8_t>& v, unsigned i) const { return value(v, aig_->outputs[i].lit); }
    uint8_t next(const std::vector<uint8_t>& v, unsigned i) const { return value(v, aig_->latches[i].next); }
    const std::vector<uint64_t>& outputSupport(unsigned i) const { return support_[aig_->outputs[i].lit >> 1]; }
    const std::vector<uint64_t>& nextSupport(unsigned i) const { return support_[aig_->latches[i].next >> 1]; }

private:
    static uint8_t value(const std::vector<uint8_t>& v, unsigned lit)
    {
        const uint8_t x = v[lit >> 1];
        return (lit & 1) ? not3(x) : x;
    }

    aiger* aig_;
    unsigned words_ = 0;
    std::vector<int> in_ap_, out_ap_;       // AIGER input/output index -> AP
    std::vector<int> ap_input_, ap_output_; // AP -> AIGER input/output index, or -1
    std::vector<std::vector<uint64_t>> support_;
};

/**
 * The product graph: a vertex per reachable product state (latch valuation,
 * automaton state) and a vertex per product edge carrying its priority, so
 * that the graph has no self-loops and every cycle passes through priorities.
 */
class Product final {
public:
    Product(const Spec& spec, const Controller& ctrl) : spec_(spec), ctrl_(ctrl) {}

    void explore()
    {
        std::deque<int> todo{productVertex(ctrl_.initialLatches(), spec_.start())};
        while (!todo.empty()) {
            const int n = todo.front();
            todo.pop_front();
            if (expanded_[n]) continue;
            expanded_[n] = true;
            const std::vector<uint8_t> latches = states_[n].first;   // copy: states_ grows below
            const int q = states_[n].second;
            std::vector<uint8_t> input(ctrl_.inputs(), VX);
            split(n, latches, q, input, todo);
        }
    }

    // FAIL iff some reachable cycle has a losing dominant priority: in every
    // SCC the dominant priority must be good, and then the SCCs that remain
    // after removing the vertices with that priority must be checked again.
    bool check() const
    {
        const int N = static_cast<int>(succ_.size());
        std::vector<std::vector<int>> work(1);
        for (int v = 0; v < N; v++) work[0].push_back(v);
        std::vector<int> region(N, 0), index(N), low(N);
        std::vector<bool> on_stack(N, false);
        int stamp = 0;
        while (!work.empty()) {
            const std::vector<int> part = std::move(work.back());
            work.pop_back();
            stamp++;
            for (int v : part) { region[v] = stamp; index[v] = -1; }
            int counter = 0;
            std::vector<int> stack;
            std::vector<std::pair<int, size_t>> call;   // iterative Tarjan: (vertex, next successor)
            for (int root : part) {
                if (index[root] != -1) continue;
                index[root] = low[root] = counter++;
                stack.push_back(root);
                on_stack[root] = true;
                call.emplace_back(root, 0);
                while (!call.empty()) {
                    const int v = call.back().first;
                    if (call.back().second < succ_[v].size()) {
                        const int w = succ_[v][call.back().second++];
                        if (region[w] != stamp) continue;
                        if (index[w] == -1) {
                            index[w] = low[w] = counter++;
                            stack.push_back(w);
                            on_stack[w] = true;
                            call.emplace_back(w, 0);
                        } else if (on_stack[w]) {
                            low[v] = std::min(low[v], index[w]);
                        }
                        continue;
                    }
                    call.pop_back();
                    if (!call.empty()) low[call.back().first] = std::min(low[call.back().first], low[v]);
                    if (low[v] != index[v]) continue;
                    std::vector<int> comp;
                    int w;
                    do {
                        w = stack.back();
                        stack.pop_back();
                        on_stack[w] = false;
                        comp.push_back(w);
                    } while (w != v);
                    if (comp.size() == 1) continue;   // trivial SCC: no self-loops in this graph
                    int dom = -1;
                    for (int y : comp) if (prio_[y] >= 0 && (dom < 0 || dominates(prio_[y], dom))) dom = prio_[y];
                    if (dom < 0) continue;
                    if (!spec_.isGood(dom)) {
                        for (int y : comp) {
                            if (prio_[y] != dom) continue;
                            printf("FAIL: a reachable cycle has dominant priority %d, e.g. through automaton state %d with latches %s\n",
                                   dom, states_[owner_[y]].second, latchString(states_[owner_[y]].first).c_str());
                            break;
                        }
                        return false;
                    }
                    std::vector<int> rest;
                    for (int y : comp) if (prio_[y] != dom) rest.push_back(y);
                    if (!rest.empty()) work.push_back(std::move(rest));
                }
            }
        }
        printf("PASS (%zu product states, %zu transitions)\n", ids_.size(), succ_.size() - ids_.size());
        return true;
    }

private:
    bool dominates(int a, int b) const { return spec_.isMax() ? a > b : a < b; }

    static std::string latchString(const std::vector<uint8_t>& latches)
    {
        std::string s;
        for (uint8_t b : latches) s += b ? '1' : '0';
        return s.empty() ? "(none)" : s;
    }

    int productVertex(const std::vector<uint8_t>& latches, int q)
    {
        const auto key = std::make_pair(latches, q);
        auto it = ids_.find(key);
        if (it != ids_.end()) return it->second;
        const int n = static_cast<int>(succ_.size());
        ids_[key] = n;
        states_.push_back(key);
        expanded_.push_back(false);
        succ_.emplace_back();
        prio_.push_back(-1);
        owner_.push_back(n);
        return n;
    }

    // Enumerate the input cubes from product state n = (latches, q): with the
    // partial input assignment, simulate and evaluate the labels in three-valued
    // logic; if the edge and the next latches are determined, add the product
    // edge, else split on an input that some undetermined signal depends on.
    void split(int n, const std::vector<uint8_t>& latches, int q, std::vector<uint8_t>& input, std::deque<int>& todo)
    {
        const std::vector<uint8_t> v = ctrl_.simulate(input, latches);
        std::vector<uint8_t> ap(spec_.data().noAPs, VX);
        for (unsigned i = 0; i < ctrl_.inputs(); i++) ap[ctrl_.inputAP(i)] = input[i];
        for (unsigned i = 0; i < ctrl_.outputs(); i++) ap[ctrl_.outputAP(i)] = ctrl_.output(v, i);
        std::vector<uint8_t> next(ctrl_.latches());
        for (unsigned i = 0; i < ctrl_.latches(); i++) next[i] = ctrl_.next(v, i);
        const int edges = spec_.edges(q);
        std::vector<uint8_t> taken(edges);
        for (int j = 0; j < edges; j++) taken[j] = spec_.eval(spec_.label(q, j), ap);

        std::vector<uint64_t> relevant(ctrl_.words(), 0);
        bool determined = true;
        for (unsigned i = 0; i < ctrl_.latches(); i++) {
            if (next[i] != VX) continue;
            determined = false;
            const std::vector<uint64_t>& s = ctrl_.nextSupport(i);
            for (unsigned w = 0; w < ctrl_.words(); w++) relevant[w] |= s[w];
        }
        for (int j = 0; j < edges; j++) {
            if (taken[j] != VX) continue;
            determined = false;
            for (int a : spec_.aps(q, j)) {
                if (ap[a] != VX) continue;
                if (ctrl_.apInput(a) >= 0) {
                    const unsigned i = ctrl_.apInput(a);
                    relevant[i / 64] |= uint64_t(1) << (i % 64);
                } else {
                    const std::vector<uint64_t>& s = ctrl_.outputSupport(ctrl_.apOutput(a));
                    for (unsigned w = 0; w < ctrl_.words(); w++) relevant[w] |= s[w];
                }
            }
        }

        if (determined) {
            // "deterministic" and "complete" are only declared in the file:
            // check them on every valuation the controller reaches.
            int edge = -1;
            for (int j = 0; j < edges; j++) {
                if (taken[j] != 1) continue;
                if (edge != -1) fail("automaton state %d has more than one edge for %s: not deterministic", q, spec_.valuationString(ap).c_str());
                edge = j;
            }
            if (edge == -1) {
                printf("FAIL: the controller reaches automaton state %d with %s, which has no edge (the automaton is not complete)\n",
                       q, spec_.valuationString(ap).c_str());
                exit(1);
            }
            const int e = static_cast<int>(succ_.size());
            succ_.emplace_back();
            prio_.push_back(spec_.prio(q, edge));
            owner_.push_back(n);
            states_.emplace_back();
            expanded_.push_back(true);
            succ_[n].push_back(e);
            const int m = productVertex(next, spec_.target(q, edge));
            succ_[e].push_back(m);
            if (!expanded_[m]) todo.push_back(m);
            return;
        }

        // Split on the first unassigned input that an undetermined signal
        // depends on (there is one: a signal is X only if some input in its
        // support is X).
        int pick = -1;
        for (unsigned i = 0; i < ctrl_.inputs() && pick < 0; i++) {
            if (input[i] == VX && ((relevant[i / 64] >> (i % 64)) & 1)) pick = static_cast<int>(i);
        }
        if (pick < 0) fail("internal error: no input to split on");
        input[pick] = 0;
        split(n, latches, q, input, todo);
        input[pick] = 1;
        split(n, latches, q, input, todo);
        input[pick] = VX;
    }

    const Spec& spec_;
    const Controller& ctrl_;
    std::map<std::pair<std::vector<uint8_t>, int>, int> ids_;
    std::vector<std::pair<std::vector<uint8_t>, int>> states_;   // by vertex id; empty for edge vertices
    std::vector<bool> expanded_;                                 // by vertex id; true for edge vertices
    std::vector<std::vector<int>> succ_;                         // all vertices
    std::vector<int> prio_;                                      // -1 for product states
    std::vector<int> owner_;                                     // vertex -> its product state (source, for edges)
};

}  // namespace

namespace {

struct Arguments {
    char** argv;
    int result;
};

void* check(void* p)
{
    Arguments& args = *static_cast<Arguments*>(p);
    char** argv = args.argv;
    FILE* f = fopen(argv[1], "r");
    if (f == nullptr) fail("cannot open %s", argv[1]);
    HoaData data;
    defaultsHoa(&data);
    const int ret = parseHoa(f, &data);
    fclose(f);
    if (ret != 0) exit(2);

    const Spec spec(data);   // may rewrite the acceptance condition of data
    const Controller ctrl(argv[2], data);
    Product product(spec, ctrl);
    product.explore();
    const bool ok = product.check();
    resetHoa(&data);
    args.result = ok ? 0 : 1;
    return nullptr;
}

}  // namespace

int main(int argc, char** argv)
{
    if (argc != 3) {
        fprintf(stderr, "usage: aigcheck automaton.ehoa controller.aag\n");
        return 2;
    }
    // The parser recurses as deep as the input is nested: give it the stack
    // that knor's Lace worker has (see main in src/knor.cpp).
    Arguments args{argv, 2};
    pthread_attr_t attr;
    pthread_t thread;
    if (pthread_attr_init(&attr) != 0 || pthread_attr_setstacksize(&attr, size_t(256) << 20) != 0
        || pthread_create(&thread, &attr, check, &args) != 0 || pthread_join(thread, nullptr) != 0) {
        fail("cannot start the checker thread");
    }
    pthread_attr_destroy(&attr);
    return args.result;
}
