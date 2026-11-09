/**
 * Copyright Tom van Dijk
 */

#ifndef KNOR_ENCODINGSEARCH_HPP
#define KNOR_ENCODINGSEARCH_HPP

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <aigcircuit.hpp>
#include <config.hpp>
#include <mealy.hpp>
#include <runstats.hpp>
#include <symgame.hpp>

extern "C" {
    #include "simplehoa.h"
}

/**
 * The encoding of the controller as an AIG: the plain encoding (as the options
 * choose it), and with don't cares also the candidates of the controller and
 * its variants, in binary and one-hot with several settings, keeping the
 * smallest; then the best binary candidates again with a variable order chosen
 * by sifting, as long as the reordering work lasts. Every candidate is
 * recorded in the statistics.
 */
class EncodingSearch final {
public:
    EncodingSearch(HoaData& data, const Config& config, RunStats& stats);

    // The smallest circuit; the variants are only used with don't cares.
    std::unique_ptr<AIGCircuit> encode(SymGame& controller, const std::vector<MealyMinimizer::Variant>& variants);

private:
    // how a game is encoded: binary or one-hot, which signals are chosen first, and which
    // simplifications of each signal are considered (see AIGEncoder::setDontCares)
    struct Setting {
        bool onehot, statesFirst, allCandidates, isopCandidate;
    };
    struct Binary {
        long size;
        size_t game;
        Setting setting;
    };

    std::vector<Setting> settings() const;
    long size(const std::unique_ptr<AIGCircuit>& circuit) const;
    // Encodes game g with the setting and keeps it if it is smaller; false if one-hot does not fit.
    bool encodeCandidate(size_t g, const Setting& s, long reorder);

    HoaData& data;
    const Config& config;
    RunStats& stats;
    std::vector<std::pair<SymGame*, std::string>> games;
    std::unique_ptr<AIGCircuit> best;
    std::vector<Binary> binaries;   // the binary candidates, to reorder the best ones
    long reorderLeft = 0;
};

#endif
