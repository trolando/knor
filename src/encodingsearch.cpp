/**
 * Copyright Tom van Dijk
 */

#include <algorithm>
#include <chrono>
#include <iostream>

#include <aigencoder.hpp>
#include <encodingsearch.hpp>

namespace {

double seconds(std::chrono::steady_clock::time_point since)
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - since).count();
}

}  // namespace

EncodingSearch::EncodingSearch(HoaData& data, const Config& config, RunStats& stats)
    : data(data), config(config), stats(stats)
{
}

long EncodingSearch::size(const std::unique_ptr<AIGCircuit>& circuit) const
{
    return circuit->getNumAnds() + (config.gatesOnly ? 0 : circuit->getNumLatches());
}

// By default only the candidates that ever gave the smallest circuit on the examples (measured
// with --dc-all --stats): all simplifications of each signal, binary also without the irredundant
// sum of products. With --dc-all every candidate, for experiments; the ablation options leave out
// some.
std::vector<EncodingSearch::Setting> EncodingSearch::settings() const
{
    const bool all = config.everyCandidate;
    const bool no_isop = !config.isopCandidate;
    const bool no_all = !config.boundCandidates;
    std::vector<Setting> result;
    for (const bool onehot : {false, true}) {
        if (!(onehot ? config.onehotCandidates : config.binaryCandidates)) continue;
        for (const bool isopCandidate : {true, false}) {
            if (no_isop && isopCandidate) continue;
            if (!all && !no_isop && onehot && !isopCandidate) continue;
            for (const bool allCandidates : {false, true}) {
                if (no_all ? allCandidates : !all && !allCandidates) continue;
                for (const bool statesFirst : {false, true}) {
                    result.push_back(Setting{onehot, statesFirst, allCandidates, isopCandidate});
                }
            }
        }
    }
    return result;
}

bool EncodingSearch::encodeCandidate(size_t g, const Setting& s, long reorder)
{
    const auto start = std::chrono::steady_clock::now();
    AIGEncoder other(data, *games[g].game);
    other.setRules(config.rules);
    if (config.isop) other.setIsop();
    if (s.onehot) other.setOneHot();
    auto candidate = other.setDontCares(s.statesFirst, s.allCandidates, s.isopCandidate).setReorder(reorder).encode();
    if (candidate == nullptr) return false;   // too many states for one-hot
    reorderLeft -= other.getReorderWork();
    RunStats::Candidate record;
    record.game = games[g].name;
    record.kind = games[g].kind;
    record.chosen = games[g].chosen;
    record.encoding = s.onehot ? "one-hot" : "binary";
    record.statesFirst = s.statesFirst;
    record.allCandidates = s.allCandidates;
    record.isopCandidate = s.isopCandidate;
    record.reorderWork = reorder;
    record.reorderUsed = other.getReorderWork();
    record.ands = candidate->getNumAnds();
    record.latches = candidate->getNumLatches();
    record.seconds = seconds(start);
    stats.candidates.push_back(std::move(record));
    if (config.verbose) {
        std::cerr << "don't cares (" << games[g].name << ", " << (s.onehot ? "one-hot, " : "binary, ")
                  << (s.statesFirst ? "state" : "outputs") << " first"
                  << (s.allCandidates ? ", all candidates" : "") << (s.isopCandidate ? "" : ", no isop")
                  << (reorder > 0 ? ", reordered" : "") << "): " << candidate->getNumAnds() << " gates and "
                  << candidate->getNumLatches() << " latches" << std::endl;
    }
    if (!s.onehot && reorder == 0) binaries.push_back(Binary{size(candidate), g, s});
    if (size(candidate) < size(best)) {
        best = std::move(candidate);
        stats.selected = static_cast<long>(stats.candidates.size()) - 1;
    }
    return true;
}

std::unique_ptr<AIGCircuit> EncodingSearch::encode(SymGame& controller,
                                                   const std::vector<MealyMinimizer::Variant>& variants)
{
    const auto start = std::chrono::steady_clock::now();
    games.clear();
    binaries.clear();
    best.reset();
    AIGEncoder encoder(data, controller);
    encoder.setRules(config.rules);
    if (config.verbose) encoder.setVerbose();
    if (config.isop) encoder.setIsop();
    if (config.onehot) encoder.setOneHot();
    if (config.sop) encoder.setSop();
    best = encoder.encode();
    RunStats::Candidate plain;
    plain.game = "controller";
    plain.kind = "controller";
    plain.encoding = "plain";
    plain.ands = best->getNumAnds();
    plain.latches = best->getNumLatches();
    plain.seconds = seconds(start);
    stats.candidates.push_back(std::move(plain));
    stats.selected = static_cast<long>(stats.candidates.size()) - 1;
    if (!config.dontCares) return std::move(best);

    // the encoding, the numbering of the states, the order of the bits and the simplifications
    // matter, in no predictable way: keep the smallest circuit (gates and latches, or with
    // --gates-only gates), also compared to the one without don't cares
    if (config.verbose) {
        std::cerr << "without don't cares: " << best->getNumAnds() << " gates and " << best->getNumLatches()
                  << " latches" << std::endl;
    }
    // not the deterministic monitor-coded games whose relaxed versions exist (those suffice), unless --dc-all
    games = {{&controller, "controller", "controller", true}};
    for (const auto& variant : variants) {
        if (config.everyCandidate || !variant.covered) {
            games.push_back({variant.game.get(), variant.name, MealyMinimizer::name(variant.kind), variant.chosen});
        }
    }
    const std::vector<Setting> all = settings();
    reorderLeft = config.reorderWork;
    for (size_t g = 0; g < games.size(); g++) {
        bool onehotFits = true;
        for (const Setting& s : all) {
            if (s.onehot && !onehotFits) continue;
            if (!encodeCandidate(g, s, 0)) onehotFits = false;
        }
    }
    // the variable order matters a lot, but sifting costs time: the best binary candidates first
    std::sort(binaries.begin(), binaries.end(), [](const Binary& a, const Binary& b) { return a.size < b.size; });
    for (size_t k = 0; k < binaries.size() && reorderLeft > 0 && !config.isop; k++) {
        encodeCandidate(binaries[k].game, binaries[k].setting,
                        std::min(reorderLeft, std::max(1L, config.reorderWork / 3)));
    }
    return std::move(best);
}
