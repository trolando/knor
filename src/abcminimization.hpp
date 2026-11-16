/**
 * Copyright Tom van Dijk
 */

#ifndef KNOR_ABCMINIMIZATION_HPP
#define KNOR_ABCMINIMIZATION_HPP

#include <chrono>
#include <vector>
#include <string>
#include <aigcircuit.hpp>
#include <memory>

extern "C" {
#include "misc/util/abc_namespaces.h"
#include "base/main/abcapis.h"
#include "base/abc/abc.h"
#include "base/main/main.h"
}

class ABCMinimization final {
public:
    explicit ABCMinimization(AIGCircuit& circuit, bool verbose = false)
        : circuit(circuit), verbose(verbose) {}
    // Rounds of drw, balance, drf and dc2 while they gain at least 1%, within the time limit (seconds; none
    // if it is 0)
    void drewrite(int seconds);
    // Rounds of ABC's compress commands (and with rrr, &rrr) while they gain (with rrr, at least 2%), within the
    // time limit (seconds; none if it is 0; &rrr gets the time that is left)
    void compress(int seconds, bool rrr = true);

private:
    AIGCircuit& circuit;
    bool verbose = false;

    static const std::vector<std::string> compressCommands;
    void executeAbcCommand(Abc_Frame_t* pAbc, const std::string& command) const;
    // the compress commands, as long as the time lasts
    void executeCompressCommands(Abc_Frame_t* pAbc, std::chrono::steady_clock::time_point end) const;
    static int getAbcNetworkSize(Abc_Frame_t* pAbc) ;
    static int getTmpFile(std::string& filename);
    void writeToAbc(Abc_Frame_t* pAbc) const;
    void readFromAbc(Abc_Frame_t* pAbc);
};

#endif //KNOR_ABCMINIMIZATION_HPP
