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

#ifndef KNOR_AIGENCODER_HPP
#define KNOR_AIGENCODER_HPP

class AIGEncoder final {
public:
    AIGEncoder(HoaData &data, SymGame &game);

    AIGEncoder& setSop() {
        this->sop = true;
        return *this;
    }

    AIGEncoder& setIsop() {
        this->isop = true;
        return *this;
    }

    AIGEncoder& setOneHot() {
        this->onehot = true;
        return *this;
    }

    AIGEncoder& setVerbose() {
        this->verbose = true;
        return *this;
    }

    std::unique_ptr<AIGCircuit> encode();

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
    unsigned int bddToAigIsop(sylvan::MTBDD bdd);
    unsigned int bddToAigCover(sylvan::ZDD bdd);       // use recursive encoding of ZDD cover (~shannon expansion)
    unsigned int bddToAigCoverSop(sylvan::ZDD cover); // use SOP encoding ("two level logic")
    void reduce(std::vector<std::vector<unsigned int>>& system, bool is_or);
    void processSOP();
    void processOnehot();
    void processBinary();
};

#endif //KNOR_AIGENCODER_HPP