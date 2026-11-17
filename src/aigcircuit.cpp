/**
 * Copyright Tom van Dijk
 */

#include <aigcircuit.hpp>
#include <stdexcept>
#include <vector>

AIGCircuit::AIGCircuit() {
    a = aiger_init();
    // we start with literal 2 (because 0/1 is reserved)
    lit = 2;
}

AIGCircuit::~AIGCircuit() {
    aiger_reset(a);
}

unsigned int AIGCircuit::makeAnd(unsigned int rhs0, unsigned int rhs1) {
    if (rhs1 < rhs0) {
        std::swap(rhs0, rhs1);
    }

    if (rhs0 == 0) return 0;
    if (rhs0 == 1) return rhs1;

    uint64_t cache_key = rhs1;
    cache_key <<= 32;
    cache_key |= rhs0;
    auto c = cache.find(cache_key);
    if (c != cache.end()) {
        return c->second;
    } else {
        aiger_add_and(a, lit, rhs0, rhs1);
        cache[cache_key] = lit;
        lit += 2;
        return lit - 2;
    }
}

void AIGCircuit::simplifyAnd(std::deque<unsigned int> &gates) {
    // TODO make non recursive
    // TODO: this should check all pairs first, before actually creating a gate??
    // for each pair of gates in gates, check the cache
    for (auto first = gates.begin(); first != gates.end(); ++first) {
        for (auto second = first + 1; second != gates.end(); ++second) {
            auto left = *first;
            auto right = *second;
            if (left > right) std::swap(left, right);
            uint64_t cache_key = right;
            cache_key <<= 32;
            cache_key |= left;
            auto c = cache.find(cache_key);
            if (c != cache.end()) {
                //gates.erase(std::remove_if(gates.begin(), gates.end(), [=](int x){return x==left or x==right;}),
                //        gates.end());
                gates.erase(second);
                gates.erase(first);
                gates.push_back(c->second);
                simplifyAnd(gates);
                return;
            }
        }
    }
}

void AIGCircuit::simplifyOr(std::deque<unsigned int> &gates) {
    // TODO make non recursive
    // TODO: this should check all pairs first, before actually creating a gate??
    // for each pair of gates in gates, check the cache
    for (auto first = gates.begin(); first != gates.end(); ++first) {
        for (auto second = first + 1; second != gates.end(); ++second) {
            int left = aiger_not(*first);
            int right = aiger_not(*second);
            if (left > right) std::swap(left, right);
            uint64_t cache_key = right;
            cache_key <<= 32;
            cache_key |= left;
            auto c = cache.find(cache_key);
            if (c != cache.end()) {
                gates.erase(second);
                gates.erase(first);
                gates.push_back(aiger_not(c->second));
                simplifyOr(gates);
                return;
            }
        }
    }
}

[[maybe_unused]] unsigned int AIGCircuit::makeMultiAnd(std::deque<unsigned int> &gates) {
    while (!gates.empty()) {
        auto last = gates.front();
        gates.pop_front();
        if (!gates.empty()) {
            auto last2 = gates.front();
            gates.pop_front();
            auto new_gate = makeAnd(last, last2);
            gates.push_back(new_gate);
        } else {
            return last;
        }
    }
    return aiger_false;
}

unsigned int AIGCircuit::makeMultiOr(std::deque<unsigned int> &gates) {
    while (!gates.empty()) {
        auto last = gates.front();
        gates.pop_front();
        if (!gates.empty()) {
            auto last2 = gates.front();
            gates.pop_front();
            auto new_gate = aiger_not(makeAnd(aiger_not(last), aiger_not(last2)));
            gates.push_back(new_gate);
        } else {
            return last;
        }
    }
    return aiger_false;
}

int AIGCircuit::writeAscii(FILE* out) {
    return aiger_write_to_file(a, aiger_ascii_mode, out);
}

int AIGCircuit::writeBinary(FILE* out) {
    return aiger_write_to_file(a, aiger_binary_mode, out);
}

void AIGCircuit::readFile(FILE* infile) {
    // into a new circuit first: when the file cannot be read, the circuit stays as it is
    aiger* read = aiger_init();
    const char *read_result = aiger_read_from_file(read, infile);
    if (read_result != nullptr) {
        const std::string message = read_result;
        aiger_reset(read);
        throw std::runtime_error("Could not read AIGER circuit from file: " + message);
    }
    aiger_delete_comments(read);
    aiger_reset(a);
    a = read;
}

bool AIGCircuit::sweep() {
    // the variables that the outputs read, also through the latches
    std::vector<bool> read(a->maxvar + 1, false);
    std::vector<unsigned> stack;
    for (unsigned i = 0; i < a->num_outputs; i++) stack.push_back(aiger_lit2var(a->outputs[i].lit));
    unsigned latches = 0, ands = 0;
    while (!stack.empty()) {
        const unsigned var = stack.back();
        stack.pop_back();
        if (var == 0 || read[var]) continue;
        read[var] = true;
        if (const aiger_and* gate = aiger_is_and(a, aiger_var2lit(var))) {
            ands++;
            stack.push_back(aiger_lit2var(gate->rhs0));
            stack.push_back(aiger_lit2var(gate->rhs1));
        } else if (const aiger_symbol* latch = aiger_is_latch(a, aiger_var2lit(var))) {
            latches++;
            stack.push_back(aiger_lit2var(latch->next));
        }
    }
    if (latches == a->num_latches && ands == a->num_ands) return false;

    aiger* swept = aiger_init();
    for (unsigned i = 0; i < a->num_inputs; i++) aiger_add_input(swept, a->inputs[i].lit, a->inputs[i].name);
    for (unsigned i = 0; i < a->num_latches; i++) {
        const aiger_symbol& latch = a->latches[i];
        if (!read[aiger_lit2var(latch.lit)]) continue;
        aiger_add_latch(swept, latch.lit, latch.next, latch.name);
        if (latch.reset != 0) aiger_add_reset(swept, latch.lit, latch.reset);
    }
    for (unsigned i = 0; i < a->num_ands; i++) {
        const aiger_and& gate = a->ands[i];
        if (read[aiger_lit2var(gate.lhs)]) aiger_add_and(swept, gate.lhs, gate.rhs0, gate.rhs1);
    }
    for (unsigned i = 0; i < a->num_outputs; i++) aiger_add_output(swept, a->outputs[i].lit, a->outputs[i].name);
    aiger_reencode(swept);
    aiger_reset(a);
    a = swept;
    // the literals were renumbered: new gates start after them, and none of the cached gates is known anymore
    lit = aiger_var2lit(a->maxvar + 1);
    cache.clear();
    return true;
}
