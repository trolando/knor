/**
 * Copyright Tom van Dijk
 */

#include <iostream>
#include <string>
#include <vector>

#include <oink/solvers.hpp>
#include <config.hpp>
#include <cxxopts.hpp>

namespace {

// with --auto-encode: the search over all winning moves, for controllers of at most so many states
constexpr int AUTO_ANY_STATES = 16;
// and with a fifth of the work of the other searches
constexpr long AUTO_ANY_SHARE = 5;
// with --auto-encode: the permissive extension of the strategy, with this much work (keeps the
// buffers' small controllers, which need about 10 million)
constexpr long AUTO_PERMISSIVE_WORK = 20000000;

// with --dc: the work (BDD nodes visited, see AIGEncoder::setReorder) for choosing the variable
// order of the binary candidates by sifting (the best candidates first), at most a third for each
constexpr long REORDER_WORK = 3000000;

}  // namespace

Config Config::parse(int argc, char** argv)
{
    try {
        cxxopts::Options opts(argv[0], "HOA synthesis using Sylvan and Oink");
        opts.custom_help("[OPTIONS...] [FILE]");
        opts.add_options()
            ("help", "Print help")
            ("sym", "Solve the parity game using the internal symbolic solver")
            ("naive", "Use the naive splitting procedure (not recommended)")
            ("explicit", "Use the explicit splitting procedure (not recommended)")
            ("real", "Only check realizability (no synthesis)")
            ("no-paritize", "Reject acceptance conditions that require extra states for conversion to parity")
            ("permissive", "Greedily extend the strategy to a winning multi-strategy (explicit solvers; slow on large games)")
            ("permissive-work", "Work of --permissive: (vertex, priority) pairs its cycle checks visit (default 0: no limit; the extension of --auto-encode without --permissive: 20 million)", cxxopts::value<long>())
            ("multi-fpi", "Use the multi-strategy of fixpoint iteration: every move winning when a vertex is decided (explicit solvers; knor's own, not Oink's solver --fpim)")
            ("multi-fpj", "Use the multi-strategy of fixpoint iteration with justifications (explicit solvers; knor's own, not Oink's solver --fpjm)")
            ("mealy-sat", "Minimise the controller as a Mealy machine with a SAT solver")
            ("mealy-any", "With --mealy-sat or --auto-encode, also search for smaller controllers that may make any winning move, without the 16-state cap of --auto-encode (explicit solvers)")
            ("mealy-work", "Work of the SAT searches of --mealy-sat (solver propagations, clauses inspected and added; about 25 million per second), as much for each wider search (with --auto-encode, a fifth for all winning moves); counted, not timed, so the result does not depend on the load (it is reproducible on one platform)", cxxopts::value<long>()->default_value("2000000000"))
            ("mealy-any-work", "Work of the search over all winning moves (default: as --mealy-work; with --auto-encode a fifth of it)", cxxopts::value<long>())
            ("bisim-game", "Apply bisimulation minimisation to the game")
            ("bisim-sol", "Apply bisimulation minimisation to the solution")
            ("bisim", "Apply bisimulation minimisation (--bisim-game and --bisim-sol)")
            ("onehot", "Use one-hot encoding for the states (recommended)")
            ("isop", "Convert BDDs to AIG using ISOP (instead of Shannon expansion)")
            ("sop", "Encode with ISOP and onehot (SOP variant of --isop --onehot)")
            ("dc", "Encode with don't cares: unused state codes and the freedom of the strategy (or with --mealy-sat, of the minimised controller), in binary and one-hot, keeping the smallest")
            ("dc-all", "With --dc, try every encoding candidate, not only those that pay off (for experiments; slower)")
            ("compress", "Compress the generated AIG using ABC")
            ("drewrite", "Compress the generated AIG using ABCs commands drw and drf")
            ("auto-encode", "Synthesize a small controller: --mealy-sat --dc with --multi-fpi (unless --multi-fpj or --sym), --permissive with a work budget, and for a small controller a short search over all winning moves")
            ("no-solve", "Do not solve, halt after constructing the parity game")
            ("print-game", "Just print the parity game (implies no-solve)")
            ("print-witness", "Print the witness parity game")
            ("print-kiss", "Print the Mealy machine in KISS format")
            ("a,write-ascii", "Write ascii AIGER file")
            ("b,write-binary", "Write binary AIGER file")
            ("v,verbose", "Be verbose")
            ("stats", "Write statistics of the run (times, minimisation, encoding candidates) as JSON to the file", cxxopts::value<std::string>())
            ;
        // for experiments: each switches off one part of --mealy-sat, --dc or --auto-encode
        opts.add_options("Ablation")
            ("no-relaxed-moves", "Encode the minimised controller's own moves only, not all its correct moves")
            ("no-renumber", "Relaxed moves keep the controller's state numbering (not breadth- and depth-first)")
            ("no-output-codes", "No state codes with output bits")
            ("no-monitor-codes", "No state codes with monitor flags")
            ("no-cube-codes", "No state codes as cubes (codes with a free bit)")
            ("no-adjacency-codes", "No state codes that give related states codes close together")
            ("no-fewer-first", "After an undecided SAT call of the minimiser, do not try fewer states first")
            ("no-unit-steps", "Over all winning moves, search the number of states with growing steps as over the strategies, not by one up to the first undecided number")
            ("no-strong-bound", "Lower bound of the minimiser only from states forced along single successors, with the initial one")
            ("no-winning-bound", "Lower bound over all winning moves only from the states that the environment forces, not also those that every winning controller visits")
            ("no-layered", "With --permissive, minimise over the extended strategy directly")
            ("no-extension", "--auto-encode without the permissive extension")
            ("no-unclarified", "No candidate with the freedom of the unclarified strategy (--dc without --mealy-sat)")
            ("no-binary-dc", "No binary candidates with don't cares")
            ("no-onehot-dc", "No one-hot candidates")
            ("no-isop-candidate", "No ISOP candidate for a signal")
            ("no-all-candidates", "For a signal only its restriction and ISOP, not its upper bound and constrain")
            ("reorder-work", "Sifting work for the variable order of the binary candidates (0: no reordering)",
             cxxopts::value<long>()->default_value(std::to_string(REORDER_WORK)))
            ("no-latch-ite", "No next-state functions as ITE on an output")
            ("no-rechoose", "Do not choose the functions again in the order found by sifting")
            ("no-marginal-choice", "Choose each function by its BDD nodes, not by the gates it adds to the circuit")
            ("no-dead-latches", "Keep the latches of state bits that stay 0")
            ("no-guard-complement", "One-hot: no \"none of the other states\" guards")
            ("no-or-gate", "Two gates for x OR f in the Shannon expansion, not one (any encoding)")
            ("gates-only", "Keep the candidate with the fewest gates, not the fewest gates and latches")
            ;
        opts.add_options("Explicit solvers")
            ("solvers", "List available solvers")
            ;

        // Oink's solvers, except any whose names knor's own options have
        std::vector<std::string> solvers;
        for (const auto & id : pg::Solvers::getSolverIDs()) {
            try {
                opts.add_options()(id, pg::Solvers::desc(id));
                solvers.push_back(id);
            } catch (const cxxopts::exceptions::option_already_exists&) {
            }
        }

        // Parse command line
        auto options = opts.parse(argc, argv);

        if (options.count("help")) {
            std::cout << opts.help() << std::endl;
            exit(0);
        }

        if (options.count("solvers")) {
            pg::Solvers::list(std::cout);
            exit(0);
        }

        if (options.count("mealy-any") && ((!options.count("mealy-sat") && !options.count("auto-encode")) || options.count("sym"))) {
            std::cerr << "error parsing options: --mealy-any requires --mealy-sat or --auto-encode, and an explicit solver (not --sym)" << std::endl;
            exit(255);
        }

        if (options.count("permissive") && options.count("sym")) {
            std::cerr << "error parsing options: --permissive requires an explicit solver (not --sym)" << std::endl;
            exit(255);
        }
        if ((options.count("multi-fpi") || options.count("multi-fpj")) && options.count("sym")) {
            std::cerr << "error parsing options: --multi-fpi and --multi-fpj require an explicit solver (not --sym)" << std::endl;
            exit(255);
        }
        if (options.count("multi-fpi") && options.count("multi-fpj")) {
            std::cerr << "error parsing options: --multi-fpi and --multi-fpj exclude each other" << std::endl;
            exit(255);
        }

        if ((options.count("dc") || options.count("auto-encode")) && (options.count("onehot") || options.count("sop"))) {
            std::cerr << "error parsing options: --dc and --auto-encode choose the encoding (not --onehot or --sop)" << std::endl;
            exit(255);
        }
        const bool dc = options.count("dc") || options.count("auto-encode");
        const bool mealy = options.count("mealy-sat") || options.count("auto-encode");
        for (const char* name : {"no-binary-dc", "no-onehot-dc", "no-isop-candidate", "no-all-candidates", "reorder-work",
                                 "no-latch-ite", "no-rechoose", "no-marginal-choice", "no-dead-latches", "no-guard-complement", "gates-only",
                                 "no-unclarified"}) {
            if (options.count(name) && !dc) {
                std::cerr << "error parsing options: --" << name << " requires --dc or --auto-encode" << std::endl;
                exit(255);
            }
        }
        for (const char* name : {"no-relaxed-moves", "no-renumber", "no-output-codes", "no-monitor-codes", "no-cube-codes",
                                 "no-adjacency-codes"}) {
            if (options.count(name) && !(dc && mealy)) {
                std::cerr << "error parsing options: --" << name << " requires --mealy-sat --dc or --auto-encode" << std::endl;
                exit(255);
            }
        }
        for (const char* name : {"no-fewer-first", "no-unit-steps", "no-strong-bound", "no-winning-bound"}) {
            if (options.count(name) && !mealy) {
                std::cerr << "error parsing options: --" << name << " requires --mealy-sat or --auto-encode" << std::endl;
                exit(255);
            }
        }
        if (options.count("mealy-any-work") && !options.count("auto-encode") && !options.count("mealy-any")) {
            std::cerr << "error parsing options: --mealy-any-work requires --auto-encode or --mealy-any" << std::endl;
            exit(255);
        }
        if (options.count("no-extension") && !options.count("auto-encode")) {
            std::cerr << "error parsing options: --no-extension requires --auto-encode" << std::endl;
            exit(255);
        }
        if (options.count("no-extension") && options.count("permissive")) {
            std::cerr << "error parsing options: --no-extension and --permissive exclude each other" << std::endl;
            exit(255);
        }
        for (const char* name : {"permissive-work", "mealy-work", "mealy-any-work", "reorder-work"}) {
            if (options.count(name) && options[name].as<long>() < 0) {
                std::cerr << "error parsing options: --" << name << " must not be negative" << std::endl;
                exit(255);
            }
        }
        if (options.count("no-layered") && !((options.count("permissive") && mealy) || options.count("auto-encode"))) {
            std::cerr << "error parsing options: --no-layered requires --permissive with --mealy-sat, or --auto-encode" << std::endl;
            exit(255);
        }
        if (options.count("auto-encode") && options.count("sym")) {
            std::cerr << "warning: --auto-encode works best with an explicit solver and --multi-fpi (the default) or --multi-fpj;"
                         " with --sym, the search over all winning moves is skipped" << std::endl;
        }

        Config config;
        if (!options.unmatched().empty()) config.input = options.unmatched()[0];
        config.verbose = options.count("verbose") > 0;
        if (options.count("stats")) config.statsPath = options["stats"].as<std::string>();

        config.paritize = options.count("no-paritize") == 0;
        config.symbolic = options.count("sym") > 0;
        for (const auto & id : solvers) {
            if (options.count(id)) config.solver = id;
        }
        config.naiveSplitting = options.count("naive") > 0;
        config.explicitSplitting = options.count("explicit") > 0;
        config.printGame = options.count("print-game") > 0;
        config.noSolve = options.count("no-solve") > 0;
        config.bisimGame = options.count("bisim") > 0 || options.count("bisim-game") > 0;
        config.bisimSolution = options.count("bisim") > 0 || options.count("bisim-sol") > 0;
        config.realizability = options.count("real") > 0;

        // --auto-encode is --mealy-sat --dc with --multi-fpi (unless --multi-fpj or --sym is given), with the permissive
        // extension (with a work budget) and the search over all winning moves for small controllers
        // (explicit solvers)
        const bool auto_encode = options.count("auto-encode") > 0;
        config.fpim = options.count("multi-fpi") > 0 || (auto_encode && !options.count("multi-fpj") && !options.count("sym"));
        config.fpjm = options.count("multi-fpj") > 0;
        config.permissive = options.count("permissive") > 0
                            || (auto_encode && !options.count("no-extension") && !options.count("sym"));
        config.permissiveWork = options.count("permissive-work") ? options["permissive-work"].as<long>()
                                : options.count("permissive") ? 0 : AUTO_PERMISSIVE_WORK;
        config.layered = options.count("no-layered") == 0;
        config.anyMoves = auto_encode || options.count("mealy-any") > 0;

        config.minimise = auto_encode || options.count("mealy-sat") > 0;
        config.dontCares = auto_encode || options.count("dc") > 0;
        MealyMinimizer::Options& minimiser = config.minimiser;
        minimiser.work = options["mealy-work"].as<long>();
        minimiser.verbose = config.verbose;
        minimiser.relax = config.dontCares;
        minimiser.variants.relaxed = options.count("no-relaxed-moves") == 0;
        minimiser.variants.renumber = options.count("no-renumber") == 0;
        minimiser.variants.outputCodes = options.count("no-output-codes") == 0;
        minimiser.variants.monitorCodes = options.count("no-monitor-codes") == 0;
        minimiser.variants.cubeCodes = options.count("no-cube-codes") == 0;
        minimiser.variants.adjacencyCodes = options.count("no-adjacency-codes") == 0;
        minimiser.fewerFirst = options.count("no-fewer-first") == 0;
        minimiser.unitSteps = options.count("no-unit-steps") == 0;
        minimiser.strongBound = options.count("no-strong-bound") == 0;
        minimiser.winningBound = options.count("no-winning-bound") == 0;
        minimiser.anyWork = options.count("mealy-any-work") ? options["mealy-any-work"].as<long>()
                            : auto_encode && options.count("mealy-any") == 0 ? minimiser.work / AUTO_ANY_SHARE
                                                                              : minimiser.work;
        minimiser.anyStates = auto_encode && options.count("mealy-any") == 0 ? AUTO_ANY_STATES : 0;

        config.isop = options.count("isop") > 0;
        config.onehot = options.count("onehot") > 0;
        config.sop = options.count("sop") > 0;
        config.rules.orGate = options.count("no-or-gate") == 0;
        config.rules.deadLatches = options.count("no-dead-latches") == 0;
        config.rules.latchSplit = options.count("no-latch-ite") == 0;
        config.rules.rechoose = options.count("no-rechoose") == 0;
        config.rules.marginalChoice = options.count("no-marginal-choice") == 0;
        config.rules.guardComplement = options.count("no-guard-complement") == 0;
        config.everyCandidate = options.count("dc-all") > 0;
        config.binaryCandidates = options.count("no-binary-dc") == 0;
        config.onehotCandidates = options.count("no-onehot-dc") == 0;
        config.isopCandidate = options.count("no-isop-candidate") == 0;
        config.boundCandidates = options.count("no-all-candidates") == 0;
        config.unclarified = options.count("no-unclarified") == 0;
        config.gatesOnly = options.count("gates-only") > 0;
        config.reorderWork = options["reorder-work"].as<long>();

        config.drewrite = options.count("drewrite") > 0;
        config.compress = options.count("compress") > 0;
        config.printKiss = options.count("print-kiss") > 0;
        config.printWitness = options.count("print-witness") > 0;
        config.output = options.count("write-binary") ? Output::binary
                        : options.count("write-ascii") ? Output::ascii : Output::none;
        return config;
    } catch (const cxxopts::exceptions::exception & e) {
        // A bad argv must exit NON-zero (10/20 are verdicts, 255 is the error
        // convention): benchmark harnesses map exit codes, and exiting 0 here
        // silently turned option mistakes into fake "success" runs.
        std::cerr << "error parsing options: " << e.what() << std::endl;
        exit(255);
    }
}

std::vector<Config::Setting> Config::settings() const
{
    std::vector<Setting> result;
    std::string stage;
    auto flag = [&](const char* name, bool value) {
        result.push_back({stage, name, value ? "yes" : "no", Setting::Type::flag});
    };
    auto number = [&](const char* name, long value) {
        result.push_back({stage, name, std::to_string(value), Setting::Type::number});
    };
    auto text = [&](const char* name, const std::string& value) {
        result.push_back({stage, name, value, Setting::Type::text});
    };
    stage = "game";
    text("input", input.empty() ? "stdin" : input);
    flag("paritize", paritize);
    text("solver", symbolic ? "sym" : solver);
    flag("naive", naiveSplitting);
    flag("explicit", explicitSplitting);
    flag("bisim-game", bisimGame);
    flag("bisim-sol", bisimSolution);
    flag("real", realizability);
    stage = "strategy";
    flag("multi-fpi", fpim);
    flag("multi-fpj", fpjm);
    flag("permissive", permissive);
    number("permissive-work", permissiveWork);
    flag("layered", layered);
    flag("any-moves", anyMoves);
    stage = "minimiser";
    const MealyMinimizer::Options& m = minimiser;
    flag("mealy-sat", minimise);
    number("work", m.work);
    number("any-work", m.anyWork);
    number("any-states", m.anyStates);
    flag("fewer-first", m.fewerFirst);
    flag("unit-steps", m.unitSteps);
    flag("strong-bound", m.strongBound);
    flag("winning-bound", m.winningBound);
    flag("relaxed-moves", m.variants.relaxed);
    flag("renumber", m.variants.renumber);
    flag("output-codes", m.variants.outputCodes);
    flag("monitor-codes", m.variants.monitorCodes);
    flag("cube-codes", m.variants.cubeCodes);
    flag("adjacency-codes", m.variants.adjacencyCodes);
    stage = "encoding";
    flag("isop", isop);
    flag("onehot", onehot);
    flag("sop", sop);
    flag("or-gate", rules.orGate);
    flag("dead-latches", rules.deadLatches);
    flag("latch-ite", rules.latchSplit);
    flag("rechoose", rules.rechoose);
    flag("marginal-choice", rules.marginalChoice);
    flag("guard-complement", rules.guardComplement);
    flag("dc", dontCares);
    flag("dc-all", everyCandidate);
    flag("binary-dc", binaryCandidates);
    flag("onehot-dc", onehotCandidates);
    flag("isop-candidate", isopCandidate);
    flag("all-candidates", boundCandidates);
    flag("unclarified", unclarified);
    flag("gates-only", gatesOnly);
    number("reorder-work", reorderWork);
    flag("drewrite", drewrite);
    flag("compress", compress);
    return result;
}

void Config::print(std::ostream& out) const
{
    std::string stage;
    for (const Setting& setting : settings()) {
        if (setting.stage != stage) {
            if (!stage.empty()) out << std::endl;
            stage = setting.stage;
            out << "configuration of the " << stage << ":";
        }
        out << " " << setting.name << "=" << setting.value;
    }
    if (!stage.empty()) out << std::endl;
}
