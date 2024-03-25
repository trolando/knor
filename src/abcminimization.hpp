/**
 * Copyright Tom van Dijk
 */

#ifndef KNOR_ABCMINIMIZATION_HPP
#define KNOR_ABCMINIMIZATION_HPP

#include <vector>
#include <string>
#include <aigmaker.hpp>

extern "C" {
#include "misc/util/abc_namespaces.h"
#include "base/main/abcapis.h"
#include "base/abc/abc.h"
#include "base/main/main.h"
}

class ABCMinimization {
private:
    AIGmaker& circuit;
    bool verbose = false;
public:
    ABCMinimization(AIGmaker& circuit, bool verbose = false) : circuit(circuit), verbose(verbose) {}
    void drewrite();
    void compress();
private:
    static const std::vector<std::string> compressCommands;
    void executeAbcCommand(Abc_Frame_t* pAbc, const std::string& command) const;
    void executeCompressCommands(Abc_Frame_t* pAbc) const;
    int getAbcNetworkSize(Abc_Frame_t* pAbc) const;
    int getTmpFile(char* tmp_filename) const;
    void writeToAbc(Abc_Frame_t* pAbc) const;
    void readFromAbc(Abc_Frame_t* pAbc);
};

#endif //KNOR_ABCMINIMIZATION_HPP
