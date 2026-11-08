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
    static constexpr int MAX_ONEHOT_STATES = 4096;

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

    /**
     * Choose the outputs and successors among all moves that game.strategies
     * and game.trans allow, using the freedom (and with the binary encoding,
     * the unused state codes) as don't cares to keep the circuit small: each
     * signal is chosen first among the outputs or the state bits (statesFirst)
     * as the smallest of its simplifications (more with allCandidates, and the
     * irredundant sum of products with isopCandidate). State bits that stay 0
     * get no latch. With the one-hot encoding, encode() returns nullptr for
     * games with more than MAX_ONEHOT_STATES states.
     */
    AIGEncoder& setDontCares(bool statesFirst = false, bool allCandidates = false, bool isopCandidate = true) {
        this->dontcares = true;
        this->statesFirst = statesFirst;
        this->allCandidates = allCandidates;
        this->isopCandidate = isopCandidate;
        return *this;
    }

    /**
     * With don't cares and the binary encoding (without ISOP), choose the order of the variables
     * for the Shannon expansion by sifting, with at most about the given work (BDD nodes visited
     * by the cost function; 0: keep the order of the game). The work is counted, not the time, so
     * that the result does not depend on the machine.
     */
    AIGEncoder& setReorder(long work) {
        this->reorderWork = work;
        return *this;
    }

    /** The work that choosing the variable order took. */
    long getReorderWork() const { return reorderUsed; }

    /**
     * Rules of the encoding that can be switched off, to measure what each
     * contributes (knor's ablation options); all are on by default.
     */
    struct Rules {
        bool orGate = true;             // x OR f as one gate in the Shannon expansion, not two
        bool deadLatches = true;        // with don't cares: no latch for a state bit that stays 0
        bool latchSplit = true;         // with don't cares: a next-state function as ITE on an output, if smaller
        bool guardComplement = true;    // one-hot with don't cares: guard a function that most states use
                                        // by "none of the other states"
        bool rechoose = true;           // binary with don't cares and reordering: choose the functions again
                                        // in the new variable order, if that is smaller
    };

    AIGEncoder& setRules(const Rules& rules) {
        this->rules = rules;
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
    bool dontcares = false; // choose moves using don't cares (binary encoding)
    bool statesFirst = false; // with don't cares: choose the next-state bits before the outputs
    bool allCandidates = false; // with don't cares: consider more simplifications of each bit
    bool isopCandidate = true; // with don't cares: consider the irredundant sum of products of each bit
    long reorderWork = 0; // with don't cares: work for choosing the variable order (binary encoding)
    long reorderUsed = 0; // the work it took
    Rules rules;

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
    bool processOnehotDontCares(int maxStates);
    void chooseWithDontCares(sylvan::MTBDD relation, sylvan::MTBDD inputVars, sylvan::MTBDD* cap_bdds,
                             sylvan::MTBDD* state_bdds);
    void rechoose(sylvan::MTBDD relation, const std::vector<bool>& live, const std::vector<uint32_t>& order,
                  uint32_t first, sylvan::MTBDD* cap_bdds, sylvan::MTBDD* state_bdds);
    unsigned int latchFunction(sylvan::MTBDD h, sylvan::MTBDD* cap_bdds, const std::vector<unsigned int>& output_lits,
                               std::deque<sylvan::MTBDD>& built);
};

#endif //KNOR_AIGENCODER_HPP