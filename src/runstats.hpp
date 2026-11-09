/**
 * Copyright Tom van Dijk
 */

#ifndef KNOR_RUNSTATS_HPP
#define KNOR_RUNSTATS_HPP

#include <chrono>
#include <string>
#include <utility>
#include <vector>

#include <config.hpp>
#include <mealy.hpp>
#include <permissive.hpp>

/**
 * Statistics of one run (--stats): the resolved configuration, the times of the
 * phases, the permissive extension, the searches of the Mealy minimiser, the
 * candidate circuits of --dc (numbered by their position, with the one the
 * encoding selected) and the circuit that was written, as one JSON object.
 */
class RunStats final {
public:
    // the version of the format: 2 added the schema, the configuration, the extension, the stop
    // reasons of the searches and the kinds, reordering work and selection of the candidates
    static constexpr int SCHEMA = 2;

    struct Candidate {
        std::string game;       // "controller" or the name of the variant
        std::string kind;       // what the game is (see MealyMinimizer::Kind), "controller" for the controller
        bool chosen = true;     // of the controller that the minimiser returned
        std::string encoding;   // "plain" (the options' encoding), "binary" or "one-hot"
        bool statesFirst = false;
        bool allCandidates = false;
        bool isopCandidate = false;
        long reorderWork = 0;   // the sifting work it was given (0: not reordered)
        long reorderUsed = 0;
        long ands = 0;
        long latches = 0;
        double seconds = 0;
    };

    struct Extension {
        bool ran = false;
        long budget = 0;        // 0: no limit
        PermissiveStrategy::Stats stats;
    };

    std::vector<Config::Setting> configuration;
    std::string verdict;   // "realizable" or "unrealizable", empty if not known
    std::vector<std::pair<std::string, double>> times;
    Extension extension;
    std::vector<MealyMinimizer::Search> searches;
    std::vector<Candidate> candidates;
    long selected = -1;    // the candidate that the encoding selected, -1 if none
    long ands = -1;        // the circuit that was written (after ABC, if used), -1 if none
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
