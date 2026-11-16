/**
 * Copyright Tom van Dijk
 */

#ifndef KNOR_ABCMINIMIZATION_HPP
#define KNOR_ABCMINIMIZATION_HPP

#include <vector>
#include <string>
#include <aigcircuit.hpp>
#include <memory>

class ABCMinimization final {
public:
    explicit ABCMinimization(AIGCircuit& circuit, bool verbose = false)
        : circuit(circuit), verbose(verbose) {}
    // Rounds of drw, balance, drf and dc2 while they gain at least 1%, within the time limit (seconds; none
    // if it is 0)
    void drewrite(int seconds);
    // Rounds of ABC's compress commands while they gain, and then with rrr rounds of &rrr and the compress commands
    // while they gain at least 2%, within the time limit (seconds; none if it is 0; &rrr gets 90% of the time that
    // is left)
    void compress(int seconds, bool rrr = true);

    // knor started by drewrite and compress to run their rounds (knor WORKER ...): runs them, and returns the
    // exit code
    static constexpr const char* WORKER = "--abc-worker";
    static int worker(int argc, char* argv[]);

private:
    AIGCircuit& circuit;
    bool verbose = false;

    // the rounds (drewrite, compress or compress-rrr) on the circuit, within the time limit (seconds)
    void minimise(int seconds, const std::string& rounds);
};

#endif //KNOR_ABCMINIMIZATION_HPP
