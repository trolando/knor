/**
 * Copyright Tom van Dijk
 */

#include <aigcircuit.hpp>
#include <abcminimization.hpp>
#include <stdexcept>

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

void AIGCircuit::simplify_and(std::deque<unsigned int> &gates) {
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
                simplify_and(gates);
                return;
            }
        }
    }
}

void AIGCircuit::simplify_or(std::deque<unsigned int> &gates) {
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
                simplify_or(gates);
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
    aiger_reset(a);
    a = aiger_init();
    const char *read_result = aiger_read_from_file(a, infile);
    if (read_result != nullptr) {
        throw std::runtime_error("Could not read AIGER circuit from file: " + std::string(read_result));
    }
    aiger_delete_comments(a);
}

void AIGCircuit::drewrite(bool verbose) {
    ABCMinimization min(*this, verbose);
    min.drewrite();
}

void AIGCircuit::compress(bool verbose) {
    ABCMinimization min(*this, verbose);
    min.compress();
}
