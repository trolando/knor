/**
 * Copyright Tom van Dijk
 */

#ifndef KNOR_RUNSTATS_HPP
#define KNOR_RUNSTATS_HPP

#include <chrono>
#include <string>
#include <utility>
#include <vector>

#include <mealy.hpp>

/**
 * Statistics of one run (--stats): the times of the phases, the searches of the
 * Mealy minimiser, the candidate circuits of --dc and the circuit that was
 * written, as one JSON object.
 */
class RunStats final {
public:
    struct Candidate {
        std::string game;       // "controller" or the name of the variant
        std::string encoding;   // "plain" (the options' encoding), "binary" or "one-hot"
        bool statesFirst = false;
        bool allCandidates = false;
        bool isopCandidate = false;
        bool reordered = false;
        long ands = 0;
        long latches = 0;
        double seconds = 0;
    };

    std::string verdict;   // "realizable" or "unrealizable", empty if not known
    std::vector<std::pair<std::string, double>> times;
    std::vector<MealyMinimizer::Search> searches;
    std::vector<Candidate> candidates;
    long ands = -1;        // the circuit that was written, -1 if none
    long latches = -1;

    void time(const std::string& phase, double seconds) { times.emplace_back(phase, seconds); }

    /**
     * Write the statistics as JSON to the given file, with the total time since
     * this object was made. Returns false if the file cannot be written.
     */
    bool write(const std::string& path) const;

private:
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
};

#endif
