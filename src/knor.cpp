/**************************************************************************
 * Copyright Tom van Dijk
 *************************************************************************/

#include <cstring> // for strcmp
#include <iostream>
#include <map>
#include <set>
#include <sys/time.h> // for gettimeofday
#include <deque>

#include <oink/game.hpp>
#include <oink/oink.hpp>
#include <config.hpp>
#include <encodingsearch.hpp>
#include <symgame.hpp>
#include <bisim.hpp>
#include <complete.hpp>
#include <abcminimization.hpp>
#include <gameconstructor.hpp>
#include <mealy.hpp>
#include <permissive.hpp>
#include <runstats.hpp>

extern "C" {
    #include <sylvan.h>
    #include "simplehoa.h"
    #include "aiger.h"
}

using namespace sylvan;

namespace {

// --stats: the statistics of the run, written when knor exits (also through exit())
struct StatsFile final {
    std::string path;
    RunStats stats;
    ~StatsFile()
    {
        if (!path.empty() && !stats.write(path)) std::cerr << "error: cannot write " << path << std::endl;
    }
} stats_file;

}  // namespace

static double wctime()
{
    struct timeval time{};
    gettimeofday(&time, nullptr);
    return (double)time.tv_sec + 1E-6 * (double)time.tv_usec;
}


VOID_TASK_0(gc_start)
{
    std::cerr << "starting garbage collection..." << std::endl;
}


VOID_TASK_0(gc_end)
{
    std::cerr << "garbage collection finished." << std::endl;
}


TASK_1(int, main_task, const Config*, _config)
{
    const Config& config = *_config;

    const bool verbose = config.verbose;
    if (verbose) config.print(std::cerr);
    stats_file.path = config.statsPath;
    RunStats& run_stats = stats_file.stats;
    run_stats.configuration = config.settings();

    // First initialize the HOA data structure
    const auto t_before_parsing = wctime();
    auto data = (HoaData*)malloc(sizeof(HoaData));
    defaultsHoa(data);

    if (config.input.empty()) {
        int ret = parseHoa(stdin, data);
        if (ret != 0) return ret;
    } else {
        const std::string& filename = config.input;
        FILE* f = fopen(filename.c_str(), "r");
        if (f == nullptr) {
            std::cerr << "file not found: " << filename << std::endl;
            return -1;
        }
        int ret = parseHoa(f, data);
        fclose(f);
        if (ret != 0) return ret;
    }

    const auto t_after_parsing = wctime();
    run_stats.time("parse", t_after_parsing - t_before_parsing);
    if (verbose) {
        std::cerr << "\033[1;37mfinished parsing automaton in " << std::fixed << (t_after_parsing - t_before_parsing) << " sec.\033[m" << std::endl;
        std::cerr << "automaton has " << data->noStates << " states." << std::endl;
    }

    // First check if the automaton is a parity automaton (after rewriting
    // acceptance conditions that are parity conditions in disguise)
    bool isMaxParity = true;
    short controllerParity = 0;
    const int originalStates = data->noStates;

    // Initialize Sylvan (the completion of the automaton compares labels as BDDs)
    sylvan_set_limits(512LL << 22, 1, 14); // should be enough (2 gigabytes)
    sylvan_init_package();
    sylvan_init_mtbdd();
    sylvan_init_zdd();
    if (verbose) sylvan_gc_hook_pregc(gc_start_CALL);
    if (verbose) sylvan_gc_hook_postgc(gc_end_CALL);

    int ret = completeHoa(data);
    if (ret != 0) return ret;
    ret = paritizeHoa(data, config.paritize);
    if (ret == 104) std::cerr << "(omit --no-paritize to allow adding memory)" << std::endl;
    if (ret != 0) return ret;
    if (verbose && data->noStates != originalStates) {
        std::cerr << "after adding memory, the automaton has " << data->noStates << " states." << std::endl;
    }
    ret = isParityGFG(data, &isMaxParity, &controllerParity);
    if (ret != 0) return ret;
    bool controllerIsOdd = controllerParity != 0;

    // Check if priorities are either all on state or all on transition, check state IDs for compatibility
    bool state_priorities = false;

    bool is_bad = false;
    for (int i=0; i<data->noStates; i++) {
        if (i != data->states[i].id) {
            std::cerr << "state " << i << " has an invalid id "  << data->states[i].id << "!" << std::endl;
            is_bad = true;
        }
        if (data->states[i].noAccSig != 0) {
            if (i == 0) {
                state_priorities = true;
            } else if (!state_priorities) {
                std::cerr << "not every state has a priority!";
                return -1;
            }
        }
    }
    if (is_bad) return -1;

    if (verbose) {
        if (state_priorities) std::cerr << "priorities are on states" << std::endl;
        else std::cerr << "priorities are on transitions" << std::endl;
    }

    const bool explicit_solver = !config.symbolic;
    const bool naive_splitting = config.naiveSplitting;
    const bool explicit_splitting = config.explicitSplitting;
    const bool write_pg = config.printGame;
    const bool no_solve = config.noSolve;
    const bool bisim_game = config.bisimGame;
    const bool bisim_sol = config.bisimSolution;

    std::unique_ptr<SymGame> sym = nullptr;
    bool realizable; // not known yet
    MTBDD anyMoves = mtbdd_false;   // with --mealy-any or --auto-encode: every move into the winning region
    mtbdd_protect(&anyMoves);
    MTBDD extendedMoves = mtbdd_false;   // with --permissive and --mealy-sat: the permissive extension
    mtbdd_protect(&extendedMoves);
    int compress_timeout = 30*60*1000; // 30 minutes x 60 seconds x 1000 milliseconds

    if (explicit_solver) {
        // Remember the start vertex
        int vstart = data->start[0];
        std::map<int, MTBDD> vertex_to_bdd;

        // Construct the explicit game
        std::unique_ptr<pg::Game> game = nullptr;
        if (naive_splitting) {
            const double t_before_splitting = wctime();
            game.reset(GameConstructor::constructGameNaive(data, isMaxParity, controllerIsOdd));
            const double t_after_splitting = wctime();
            run_stats.time("construct", t_after_splitting - t_before_splitting);
            if (verbose) {
                std::cerr << "\033[1;37mfinished constructing game in " << std::fixed << (t_after_splitting - t_before_splitting) << " sec.\033[m" << std::endl;
            }
        } else if (explicit_splitting) {
            const double t_before_splitting = wctime();
            game.reset(GameConstructor::constructGame(data, isMaxParity, controllerIsOdd));
            const double t_after_splitting = wctime();
            run_stats.time("construct", t_after_splitting - t_before_splitting);
            if (verbose) {
                std::cerr << "\033[1;37mfinished constructing game in " << std::fixed << (t_after_splitting - t_before_splitting) << " sec.\033[m" << std::endl;
            }
        } else {
            const double t_1 = wctime();
            vstart = 0; // always set to 0 by constructSymGame
            sym = GameConstructor::constructSymGame(data, isMaxParity, controllerIsOdd);
            const double t_2 = wctime();
            run_stats.time("construct", t_2 - t_1);
            if (verbose) std::cerr << "\033[1;37mfinished constructing symbolic game in " << std::fixed << (t_2 - t_1) << " sec.\033[m" << std::endl;
            if (bisim_game) {
                const double t_before = wctime();
                MTBDD partition = CALL(min_lts_strong, sym.get(), false);
                mtbdd_protect(&partition);
                CALL(minimize, sym.get(), partition, verbose);
                mtbdd_unprotect(&partition);
                if (verbose) {
                     std::cerr << "after bisimulation minimisation: " << count_blocks() << " blocks." << std::endl;
                }
                const double t_after = wctime();
                run_stats.time("bisimulation of the game", t_after - t_before);
                if (verbose) std::cerr << "\033[1;37mfinished bisimulation minimisation of game in " << std::fixed << (t_after - t_before) << " sec.\033[m" << std::endl;
            }
            const double t_3 = wctime();
            game = sym->toExplicit(vertex_to_bdd);
            const double t_4 = wctime();
            run_stats.time("explicit game", t_4 - t_3);
            if (verbose) {
                std::cerr << "\033[1;37mfinished constructing explicit game in " << std::fixed << (t_4 - t_3) << " sec.\033[m" << std::endl;
            }
        }

        if (verbose) {
            std::cerr << "constructed game has " << game->vertexcount() << " vertices and " << game->edgecount() << " edges." << std::endl;
        }

        if (write_pg) {
            // in case we want to write the file to PGsolver file format...
            game->set_label(vstart, "initial");
            game->write_pgsolver(std::cout);
            exit(0);
        }

        if (no_solve) {
            exit(0);
        }

        // we sort now, so we can track the initial state
        double begin = wctime();
        std::vector<int> mapping(game->vertexcount());
        game->sort(mapping.data());

        // OK, fire up the engine
        std::stringstream log;

        const std::string& solver = config.solver;

        pg::Oink engine(*game, verbose ? std::cerr : log);
        engine.setTrace(0); //verbose ? 1 : 0); actually don´t -- maybe add a 2nd verbosity level later
        engine.setRenumber();
        engine.setSolver(solver);
        engine.setWorkers(-1);

        // and run the solver
        engine.run();
        double end = wctime();
        run_stats.time("solve", end - begin);

        // report how long it all took
        if (verbose) {
            std::cerr << "\033[1;37mfinished solving game in " << std::fixed << (end - begin) << " sec.\033[m" << std::endl;
        }

        // undo the sorting, of the game and of its solution ## FOR SOME REASON, THIS CAN BE SLOW!
        pg::permute(*game, engine.solution(), mapping.data());
        mapping.clear(); // erase result, we don't need it anymore
        const pg::Solution& solution = engine.solution();

        realizable = solution.winner(vstart) == 0;

        // finally, check if the initial vertex is won by controller or environment, and unless
        // only the verdict is asked for, get the strategy
        if (realizable && !config.realizability) {
            if (sym != nullptr && config.anyMoves) {
                std::map<MTBDD, std::vector<MTBDD>> str;  // cap to winning priostates
                for (int v=0; v<game->vertexcount(); v++) {
                    if (game->owner(v) != 0 || solution.winner(v) != 0) continue;
                    for (const int* e = game->outs(v); *e != -1; e++) {
                        if (solution.winner(*e) == 0) str[vertex_to_bdd.at(v)].push_back(vertex_to_bdd.at(*e));
                    }
                }
                if (sym->applyMultiStrategy(str)) anyMoves = sym->strategies;
            }
            const bool fpim = config.fpim;
            const bool fpjm = config.fpjm;
            if (sym != nullptr && (fpim || fpjm || config.permissive)) {
                auto allowed = PermissiveStrategy::positional(*game, solution);
                if (fpim || fpjm) {
                    const double t_before = wctime();
                    const char* name = fpim ? "fixpoint iteration" : "fixpoint iteration with justifications";
                    std::vector<int> winner;
                    auto multi = fpim ? PermissiveStrategy::fpim(*game, winner) : PermissiveStrategy::fpjm(*game, winner);
                    bool agrees = true;
                    for (int v=0; v<game->vertexcount(); v++) agrees = agrees && winner[v] == solution.winner(v);
                    if (agrees && PermissiveStrategy::verify(*game, solution, vstart, multi)) {
                        allowed = std::move(multi);
                    } else {
                        // cannot happen; keep the positional strategy rather than risk a wrong controller
                        std::cerr << "error: " << name << " disagrees with " << solver << ", using its strategy" << std::endl;
                    }
                    run_stats.time(fpim ? "fpim" : "fpjm", wctime() - t_before);
                    if (verbose) {
                        long moves = 0;
                        for (const auto& a : allowed) moves += static_cast<long>(a.size());
                        std::cerr << name << ": " << moves << " moves of the controller" << std::endl;
                        std::cerr << "\033[1;37mfinished " << name << " in " << std::fixed << (wctime() - t_before) << " sec.\033[m" << std::endl;
                    }
                }
                std::vector<std::vector<int>> extended;  // the permissive extension of allowed, if any
                if (config.permissive) {
                    const double t_before = wctime();
                    PermissiveStrategy::Stats stats;
                    extended = PermissiveStrategy::extend(*game, solution, vstart, allowed, config.permissiveWork, stats);
                    run_stats.extension = RunStats::Extension{true, config.permissiveWork, stats};
                    if (!PermissiveStrategy::verify(*game, solution, vstart, extended)) {
                        // cannot happen; keep the strategy that was extended rather than risk a wrong controller
                        std::cerr << "error: the multi-strategy is not winning, using the strategy before extending it" << std::endl;
                        extended.clear();
                    }
                    run_stats.time("extension", wctime() - t_before);
                    if (verbose) {
                        std::cerr << "multi-strategy: " << stats.added << " of " << stats.candidates
                                  << " further edges allowed at " << stats.vertices << " reachable vertices, work "
                                  << stats.work << (stats.stopped ? " (stopped)" : "") << std::endl;
                        std::cerr << "\033[1;37mfinished extending the strategy in " << std::fixed << (wctime() - t_before) << " sec.\033[m" << std::endl;
                    }
                }
                auto apply = [&](const std::vector<std::vector<int>>& moves) {
                    std::map<MTBDD, std::vector<MTBDD>> str;  // cap to allowed priostates
                    for (int v=0; v<game->vertexcount(); v++) {
                        for (int w : moves[v]) str[vertex_to_bdd.at(v)].push_back(vertex_to_bdd.at(w));
                    }
                    if (!sym->applyMultiStrategy(str)) {
                        std::cerr << "cannot apply strategy!" << std::endl;
                    }
                };
                if (!extended.empty() && config.minimise && config.layered) {
                    // the minimiser searches the extension after the strategy (see MealyMinimizer::minimize)
                    apply(extended);
                    extendedMoves = sym->strategies;
                    apply(allowed);
                } else {
                    apply(extended.empty() ? allowed : extended);
                }
            } else if (sym != nullptr) {
                // now get the strategy from Oink...
                std::map<MTBDD, MTBDD> str;  // cap to priostate
                for (int v=0; v<game->vertexcount(); v++) {
                    // good cap states are owned and won by 0
                    if (game->owner(v) == 0 && solution.winner(v) == 0) {
                        auto a = vertex_to_bdd.at(v);
                        auto b = vertex_to_bdd.at(solution.strategy(v));
                        str[a] = b;
                    }
                }

                if (!sym->applyStrategy(str)) {
                    std::cerr << "cannot apply strategy!" << std::endl;
                }
            }
        }
    } else {
        // Construct the game
        const double t_before_construct = wctime();
        sym = GameConstructor::constructSymGame(data, isMaxParity, controllerIsOdd);
        const double t_after_construct = wctime();
        run_stats.time("construct", t_after_construct - t_before_construct);

        if (verbose) {
            std::cerr << "\033[1;37mfinished constructing symbolic game in " << std::fixed << (t_after_construct - t_before_construct) << " sec.\033[m" << std::endl;
        }

        if (bisim_game) {
            const double t_before = wctime();
            MTBDD partition = CALL(min_lts_strong, sym.get(), false);
            mtbdd_protect(&partition);
            CALL(minimize, sym.get(), partition, verbose);
            mtbdd_unprotect(&partition);
            if (verbose) {
                std::cerr << "after bisimulation minimisation: " << count_blocks() << " blocks." << std::endl;
            }
            const double t_after = wctime();
            run_stats.time("bisimulation of the game", t_after - t_before);
            if (verbose) std::cerr << "\033[1;37mfinished bisimulation minimisation of game in " << std::fixed << (t_after - t_before) << " sec.\033[m" << std::endl;
        }

        if (write_pg) {
            // in case we want to write the file to PGsolver file format...
            std::map<int, MTBDD> vertex_to_bdd;
            auto pg = sym->toExplicit(vertex_to_bdd);
            pg->write_pgsolver(std::cout);
            exit(0);
        }

        if (no_solve) {
            exit(0);
        }

        const double t_before_solve = wctime();
        realizable = sym->solve(verbose);
        const double t_after_solve = wctime();
        run_stats.time("solve", t_after_solve - t_before_solve);

        if (verbose) {
            std::cerr << "\033[1;37mfinished solving game in " << std::fixed << (t_after_solve - t_before_solve) << " sec.\033[m" << std::endl;
        }
    }

    run_stats.verdict = realizable ? "realizable" : "unrealizable";

    if (config.realizability) {
        if (realizable) {
            std::cout << "REALIZABLE" << std::endl;
            if (verbose) {
                std::cerr << "\033[1;37mtotal time was " << std::fixed << (wctime() - t_before_parsing) << " sec.\033[m" << std::endl;
            }
            exit(10);
        } else {
            std::cout << "UNREALIZABLE" << std::endl;
            if (verbose) {
                std::cerr << "\033[1;37mtotal time was " << std::fixed << (wctime() - t_before_parsing) << " sec.\033[m" << std::endl;
            }
            exit(20);
        }
    }

    if (!realizable) {
        if (verbose) {
            std::cerr << "\033[1;37mtotal time was " << std::fixed << (wctime() - t_before_parsing) << " sec.\033[m"
                      << std::endl;
            std::cerr << "\033[1;31mgame is unrealizable!\033[m" << std::endl;
            sylvan_stats_report(stdout);
        }
        mtbdd_unprotect(&anyMoves);
        mtbdd_unprotect(&extendedMoves);
        return 20;
    }

    if (verbose) std::cerr << "\033[1;38;5;10mgame is realizable!\033[m" << std::endl;

    if (naive_splitting or explicit_splitting) {
        std::cerr << "--naive and --explicit are incompatible with generating the AIG!" << std::endl;
        exit(10);
    }

    if (verbose) {
        // sym->print_trans(true);
        // sym->print_strategies();
    }

    // with --dc, the minimised controller also comes with other state codes and all its correct
    // moves, for the encoder to choose from; such games are not post-processed or bisimulated
    std::vector<MealyMinimizer::Variant> variants;
    bool minimised = false;
    if (config.minimise) {
        const auto t_before = wctime();
        MealyMinimizer::Options minimise = config.minimiser;
        minimise.extended = extendedMoves;
        minimise.anyMoves = anyMoves;
        auto minimized = MealyMinimizer::minimize(*sym, minimise);
        run_stats.searches = minimized.searches;
        run_stats.time("minimisation", wctime() - t_before);
        if (minimized.game != nullptr) {
            sym = std::move(minimized.game);
            variants = std::move(minimized.variants);
            minimised = true;
            if (verbose) std::cerr << "minimised controller has " << minimized.states << " states." << std::endl;
        }
        if (verbose) {
            std::cerr << "\033[1;37mfinished mealy minimisation in " << std::fixed << (wctime() - t_before) << " sec.\033[m" << std::endl;
        }
    }
    mtbdd_unprotect(&anyMoves);
    mtbdd_unprotect(&extendedMoves);

    // without a minimised controller, --dc may choose among the moves of the strategy itself, rather
    // than the lowest output that postprocess picks: on the explicit route, the strategy is a winning
    // multi-strategy (the freedom of the symbolic solver need not be), so every choice is correct
    if (config.dontCares && !minimised && explicit_solver && config.unclarified) {
        auto relational = std::make_unique<SymGame>(sym->statebits, sym->priobits, sym->uap_count, sym->cap_count,
                                                     sym->maxprio);
        relational->trans = sym->trans;
        relational->strategies = sym->strategies;
        relational->restrictToReachable(false);
        variants.push_back({std::move(relational), "unclarified strategy", MealyMinimizer::Kind::unclarified});
    }

    {
        const auto t_before = wctime();
        sym->postprocess(verbose);
        const auto t_after = wctime();
        run_stats.time("postprocess", t_after - t_before);

        if (verbose) {
            std::cerr << "\033[1;37mfinished post processing in " << std::fixed << (t_after - t_before) << " sec.\033[m" << std::endl;
        }
    }

    if (verbose) {
        // sym->print_trans();
        // sym->print_strategies();
    }

    if (bisim_sol) {
        const auto t_before = wctime();
        MTBDD partition = CALL(min_lts_strong, sym.get(), true);
        mtbdd_protect(&partition);
        if (verbose) {
            // CALL(print_partition, sym, partition);
        }
        CALL(minimize, sym.get(), partition, verbose);
        mtbdd_unprotect(&partition);
        const auto t_after = wctime();
        run_stats.time("bisimulation of the solution", t_after - t_before);
        if (verbose) std::cerr << "\033[1;37mfinished bisimulation minimisation of solution in " << std::fixed << (t_after - t_before) << " sec.\033[m" << std::endl;
    }

    if (config.printKiss) {
        sym->print_kiss(true);
        if (verbose) {
            std::cerr << "\033[1;37mtotal time was " << std::fixed << (wctime() - t_before_parsing) << " sec.\033[m" << std::endl;
        }
        exit(10);
    }

    /**
     * maybe print witness parity game, which should be fully won by even
     */

    if (config.printWitness) {
        auto pargame = sym->strategyToPG();
        pargame->write_pgsolver(std::cout);
        if (verbose) {
            std::cerr << "\033[1;37mtotal time was " << std::fixed << (wctime() - t_before_parsing) << " sec.\033[m" << std::endl;
        }
        exit(10);
    }

    const double t_before_enc = wctime();
    auto circuit = EncodingSearch(*data, config, run_stats).encode(*sym, variants);
    const double t_after_enc = wctime();
    run_stats.time("encoding", t_after_enc - t_before_enc);
    if (verbose) std::cerr << "\033[1;37mfinished encoding in " << std::fixed << (t_after_enc - t_before_enc) << " sec.\033[m" << std::endl;

    /**
     * maybe compress with ABC
     */
    if (config.drewrite) {
        if (verbose) std::cerr << "size of AIG before drw+drf: " << circuit->getNumAnds() << " gates." << std::endl;
        const double t_before = wctime();
        ABCMinimization(*circuit, verbose).drewrite(compress_timeout);
        const double t_after = wctime();
        if (verbose) std::cerr << "size of AIG after drw+drf: " << circuit->getNumAnds() << " gates." << std::endl;
        run_stats.time("drw+drf", t_after - t_before);
        if (verbose) std::cerr << "\033[1;37mfinished drw+drf in " << std::fixed << (t_after - t_before) << " sec.\033[m" << std::endl;
    }

    if (config.compress) {
        if (verbose) std::cerr << "size of AIG before compression: " << circuit->getNumAnds() << " gates." << std::endl;
        const double t_before = wctime();
        ABCMinimization(*circuit, verbose).compress(compress_timeout);
        const double t_after = wctime();
        if (verbose) std::cerr << "size of AIG after compression: " << circuit->getNumAnds() << " gates." << std::endl;
        run_stats.time("compression", t_after - t_before);
        if (verbose) std::cerr << "\033[1;37mfinished compression in " << std::fixed << (t_after - t_before) << " sec.\033[m" << std::endl;
    }

    if (verbose) {
        std::cerr << "final size of AIG: " << circuit->getNumAnds() << " gates." << std::endl;
    }
    run_stats.ands = circuit->getNumAnds();
    run_stats.latches = circuit->getNumLatches();

    if (config.output == Config::Output::binary) {
        circuit->writeBinary(stdout);
    } else if (config.output == Config::Output::ascii) {
        circuit->writeAscii(stdout);
    }

    if (verbose) {
        std::cerr << "\033[1;37mtotal time was " << std::fixed << (wctime() - t_before_parsing) << " sec.\033[m" << std::endl;
        sylvan_stats_report(stdout);
    }

    return 10;

    // We don't need Sylvan anymore at this point
    // sylvan_quit();
    // free HOA allocated data structure
    // resetHoa(data);
}


/**
 * The main function
 */
int
main(int argc, char* argv[])
{
    const Config config = Config::parse(argc, argv);

    // Initialize Lace, only 1 worker. The worker runs everything, and its program stack must hold
    // the recursion over deeply nested acceptance conditions (parity games with up to 100000
    // priorities, see MAX_NESTING in simplehoa.cpp); only the part that is used takes memory.
    lace_set_stacksize(size_t(256) << 20);
    lace_start(1, 1024*1024*2); // initialize Lace, but sequentially
                                // also get a large enough task size... (2M tasks) for PSI!

    int res = RUN(main_task, &config);

    lace_stop();

    return res;
}
