/**
 * Copyright Tom van Dijk
 */

#include <sylvan.h>
#include <oink.hpp>
#include <deque>
#include <symgame.hpp>

extern "C" {
    #include "aiger.h"
    #include "simplehoa.h"
    #include "misc/util/abc_namespaces.h"
    #include "base/main/abcapis.h"
    #include "base/abc/abc.h"
    #include "base/main/main.h"
}

#pragma once

class AIGmaker {
private:
    aiger *a;
    HoaData *data;
    SymGame *game;

    bool isop = false; // use ISOP
    bool verbose = false;
    bool onehot = false; // use onehot encoding
    
    int lit; // current next literal

    int* uap_to_lit; // the input literal for each uncontrolled AP
    std::map<int, int> state_to_lit; // the latch literal for each state bit / onehot state
    char** caps; // labels for controlled APs
    std::map<uint32_t, int> var_to_lit; // translate BDD variable (uap/state) to AIGER literal

    std::map<sylvan::MTBDD, int> mapping; // map MTBDD to AIGER literal
    std::map<uint64_t, int> cache; // cache for ands

    int makeand(int rhs0, int rhs1);
    int bdd_to_aig(sylvan::MTBDD bdd);           // use recursive encoding of BDD (shannon expanion)
    int bdd_to_aig_isop(sylvan::MTBDD bdd);
    int bdd_to_aig_cover(sylvan::ZDD bdd);       // use recursive encoding of ZDD cover (~shannon expansion)
    int bdd_to_aig_cover_sop(sylvan::ZDD cover); // use SOP encoding ("two level logic")
    void simplify_and(std::deque<int> &gates);
    void simplify_or(std::deque<int> &gates);

    int make_and(std::deque<int> &gates);
    int make_or(std::deque<int> &gates);

    static const std::vector<std::string> compressCommands;
    void executeAbcCommand(Abc_Frame_t* pAbc, const std::string command) const;
    void executeCompressCommands(Abc_Frame_t* pAbc) const;
    int getAbcNetworkSize(Abc_Frame_t* pAbc) const;
    int getTmpFile(char* tmp_filename) const;
    void writeToAbc(Abc_Frame_t* pAbc) const;
    void readFromAbc(Abc_Frame_t* pAbc);

public:
    AIGmaker(HoaData *data, SymGame *game);
    ~AIGmaker();

    void setIsop()
    {
        this->isop = true;
    }

    void setOneHot()
    {
        this->onehot = true;
    }

    void setVerbose()
    {
        this->verbose = true;
    }

    long getNumAnds()
    {
        return this->a->num_ands;
    }

    void process();
    void writeAscii(FILE* out);
    void writeBinary(FILE* out);

    void compress();
};

