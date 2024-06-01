/**
 * Copyright Tom van Dijk
 */

#ifndef KNOR_AIGCIRCUIT_HPP
#define KNOR_AIGCIRCUIT_HPP

#include <deque>
#include <map>
#include <cstdint>
#include <string>
#include <iostream>

extern "C" {
#include <aiger.h>
}

class AIGCircuit final {
public:
    AIGCircuit();
    ~AIGCircuit();

    unsigned int makeAnd(unsigned int rhs0, unsigned int rhs1);
    [[maybe_unused]] unsigned int makeMultiAnd(std::deque<unsigned int> &gates);
    [[maybe_unused]] unsigned int makeMultiOr(std::deque<unsigned int> &gates);

    void simplify_and(std::deque<unsigned int> &gates);
    void simplify_or(std::deque<unsigned int> &gates);

    long getNumAnds() {
        return this->a->num_ands;
    }

    long getNumLatches() {
        return this->a->num_latches;
    }

    unsigned int makeInput(const std::string& label) {
        aiger_add_input(a, lit, label.c_str());
        lit += 2;
        return lit - 2;
    }

    unsigned int makeLatch() {
        lit += 2;
        return lit - 2;
    }

    void setLatch(unsigned int latch, unsigned int gate, const std::string& label) {
        aiger_add_latch(a, latch, gate, label.c_str());
    }

    void makeOutput(unsigned int gate, const std::string& label) {
        aiger_add_output(a, gate, label.c_str());
    }

    int writeAscii(FILE* outfile);
    int writeBinary(FILE* outfile);
    void readFile(FILE* infile);

    void drewrite(bool verbose);
    void compress(bool verbose);

private:
    aiger *a;
    std::map<uint64_t, unsigned int> cache; // cache for ands
    unsigned int lit; // current next literal
};

#endif //KNOR_AIGCIRCUIT_HPP
