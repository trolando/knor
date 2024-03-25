/**
 * Copyright Tom van Dijk
 */

#include <sylvan.h>
#include <oink/oink.hpp>
#include <deque>
#include <symgame.hpp>
#include <aigcircuit.hpp>
#include <sylvan_mtbdd.h>

extern "C" {
    #include "aiger.h"
    #include "simplehoa.h"
}

#pragma once

class AIGmaker final {
public:
    AIGmaker(HoaData *data, SymGame *game);

    AIGmaker& setSop()
    {
        this->sop = true;
        return *this;
    }

    AIGmaker& setIsop()
    {
        this->isop = true;
        return *this;
    }

    AIGmaker& setOneHot()
    {
        this->onehot = true;
        return *this;
    }

    AIGmaker& setVerbose()
    {
        this->verbose = true;
        return *this;
    }

    std::unique_ptr<AIGCircuit> process();
private:
    HoaData &data;
    SymGame &game;
    std::unique_ptr<AIGCircuit> circuit;

    bool sop = false; // use ISOP onehot SOP
    bool isop = false; // use ISOP
    bool verbose = false;
    bool onehot = false; // use onehot encoding

    std::vector<unsigned int> uap_to_lit; // store the input literals for each uncontrolled AP
    std::vector<std::string> cap_labels; // store the labels for controlled APs
    std::map<int, unsigned int> state_to_lit; // the latch literal for each state bit / onehot state
    std::map<uint32_t, unsigned int> bddvar_to_lit; // translate BDD variable (uap/state) to AIGER literal

    unsigned int bddToAigRecursive(sylvan::MTBDD bdd);           // use recursive encoding of BDD (shannon expanion)
    unsigned int bdd_to_aig_isop(sylvan::MTBDD bdd);
    unsigned int bddToAigCover(sylvan::ZDD bdd);       // use recursive encoding of ZDD cover (~shannon expansion)
    unsigned int bdd_to_aig_cover_sop(sylvan::ZDD cover); // use SOP encoding ("two level logic")
    void reduce(std::vector<std::vector<unsigned int>>& system, bool is_or);
    void processSOP();
    void processOnehot();
    void processBinary();
};

