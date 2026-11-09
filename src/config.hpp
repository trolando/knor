/**
 * Copyright Tom van Dijk
 */

#ifndef KNOR_CONFIG_HPP
#define KNOR_CONFIG_HPP

#include <iosfwd>
#include <string>

#include <aigencoder.hpp>
#include <mealy.hpp>

/**
 * The configuration of a run: the command line parsed, checked and resolved
 * once (what --auto-encode implies, the default budgets, the solver), so that
 * the stages see what to do rather than the options that say it.
 */
struct Config final {
    std::string input;                  // the automaton (empty: standard input)
    bool verbose = false;
    std::string statsPath;              // --stats (empty: none)

    // the parity game and its solution
    bool paritize = false;
    bool symbolic = false;              // the symbolic solver (--sym), else Oink's
    std::string solver = "tl";          // Oink's solver
    bool naiveSplitting = false;
    bool explicitSplitting = false;
    bool printGame = false;
    bool noSolve = false;
    bool bisimGame = false;
    bool bisimSolution = false;
    bool realizability = false;         // only the verdict (--real)

    // the strategy, of Oink's solution
    bool fpim = false;                  // the multi-strategy of fixpoint iteration
    bool fpjm = false;                  // the multi-strategy of fixpoint iteration with justifications
    bool permissive = false;            // the permissive extension of the strategy
    long permissiveWork = 0;            // its work (0: no limit)
    bool layered = true;                // with permissive and minimise: search the extension after the strategy
    bool anyMoves = false;              // all moves into the winning region, for the search over them

    // the minimisation as a Mealy machine (the moves of the options are filled in by the run)
    bool minimise = false;
    MealyMinimizer::Options minimiser;

    // the encoding
    bool isop = false;                  // the plain encoding with ISOP
    bool onehot = false;                // the plain encoding one-hot
    bool sop = false;                   // the plain encoding with ISOP and one-hot
    AIGEncoder::Rules rules;
    bool dontCares = false;             // also the candidates with don't cares, keeping the smallest
    bool everyCandidate = false;        // every candidate, also those that never pay off (--dc-all)
    bool binaryCandidates = true;
    bool onehotCandidates = true;
    bool isopCandidate = true;          // for a signal also its ISOP
    bool boundCandidates = true;        // for a signal also its upper bound and constrain
    bool unclarified = true;            // without a minimised controller, a candidate with the strategy's freedom
    bool gatesOnly = false;             // keep the candidate with the fewest gates, not gates and latches
    long reorderWork = 0;               // sifting work for the binary candidates

    // after the encoding
    bool drewrite = false;
    bool compress = false;
    bool printKiss = false;
    bool printWitness = false;
    enum class Output { none, ascii, binary } output = Output::none;

    // Parse and check the command line; exits on --help, --solvers and errors (255).
    static Config parse(int argc, char** argv);

    // The resolved configuration, a line per stage, named after the options.
    void print(std::ostream& out) const;
};

#endif
