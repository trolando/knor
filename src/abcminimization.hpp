/**
 * Copyright Tom van Dijk
 */

#ifndef KNOR_ABCMINIMIZATION_HPP
#define KNOR_ABCMINIMIZATION_HPP

#include <vector>
#include <string>
#include <aigcircuit.hpp>

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
    void drewrite(int timeout);
    void compress(int timeout);

private:
    AIGCircuit& circuit;
    bool verbose = false;

    static const std::vector<std::string> compressCommands;
    void executeAbcCommand(Abc_Frame_t* pAbc, const std::string& command) const;
    void executeCompressCommands(Abc_Frame_t* pAbc) const;
    static int getAbcNetworkSize(Abc_Frame_t* pAbc) ;
    static int getTmpFile(char* tmp_filename) ;
    void writeToAbc(Abc_Frame_t* pAbc) const;
    void readFromAbc(Abc_Frame_t* pAbc);
};

#endif //KNOR_ABCMINIMIZATION_HPP
