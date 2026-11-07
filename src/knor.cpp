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
#include <oink/solvers.hpp>
#include <cxxopts.hpp>
#include <symgame.hpp>
#include <bisim.hpp>
#include <aigencoder.hpp>
#include <abcminimization.hpp>
#include <gameconstructor.hpp>
#include <mealy.hpp>
#include <permissive.hpp>

extern "C" {
    #include <sylvan.h>
    #include "simplehoa.h"
    #include "aiger.h"
}

using namespace sylvan;

static double wctime()
{
    struct timeval time{};
    gettimeofday(&time, nullptr);
    return (double)time.tv_sec + 1E-6 * (double)time.tv_usec;
}


cxxopts::ParseResult
handleOptions(int &argc, char**& argv)
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
            ("paritize", "Convert any acceptance condition (e.g. generalized Buchi) to parity, adding memory")
            ("permissive", "Greedily extend the strategy to a winning multi-strategy (explicit solvers; slow on large games)")
            ("fpim", "Use the multi-strategy of fixpoint iteration: every move winning when a vertex is decided (explicit solvers)")
            ("fpjm", "Use the multi-strategy of fixpoint iteration with justifications (explicit solvers)")
            ("mealy-sat", "Minimise the controller as a Mealy machine with a SAT solver")
            ("mealy-any", "With --mealy-sat, also search for smaller controllers that may make any winning move (explicit solvers)")
            ("mealy-time", "Seconds for the SAT searches of --mealy-sat, and as many for each wider search; the rest of the minimisation is not bounded", cxxopts::value<double>()->default_value("60"))
            ("bisim-game", "Apply bisimulation minimisation to the game")
            ("bisim-sol", "Apply bisimulation minimisation to the solution")
            ("bisim", "Apply bisimulation minimisation (--bisim-game and --bisim-sol)")
            ("onehot", "Use one-hot encoding for the states (recommended)")
            ("isop", "Convert BDDs to AIG using ISOP (instead of Shannon expansion)")
            ("sop", "Encode with ISOP and onehot (SOP variant of --isop --onehot)")
            ("compress", "Compress the generated AIG using ABC")
            ("drewrite", "Compress the generated AIG using ABCs commands drw and drf")
            ("best", "Try all combinations of bisim and isop and write the smallest AIG")
            ("no-solve", "Do not solve, halt after constructing the parity game")
            ("print-game", "Just print the parity game (implies no-solve)")
            ("print-witness", "Print the witness parity game")
            ("print-kiss", "Print the Mealy machine in KISS format")
            ("a,write-ascii", "Write ascii AIGER file")
            ("b,write-binary", "Write binary AIGER file")
            ("v,verbose", "Be verbose")
            ;
        opts.add_options("Explicit solvers")
            ("solvers", "List available solvers")
            ;

        // Add solvers
        for (const auto & id : pg::Solvers::getSolverIDs()) {
            opts.add_options()(id, pg::Solvers::desc(id));
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

        if (options.count("mealy-any") && (!options.count("mealy-sat") || options.count("sym"))) {
            std::cerr << "error parsing options: --mealy-any requires --mealy-sat and an explicit solver (not --sym)" << std::endl;
            exit(255);
        }

        if (options.count("permissive") && options.count("sym")) {
            std::cerr << "error parsing options: --permissive requires an explicit solver (not --sym)" << std::endl;
            exit(255);
        }
        if ((options.count("fpim") || options.count("fpjm")) && options.count("sym")) {
            std::cerr << "error parsing options: --fpim and --fpjm require an explicit solver (not --sym)" << std::endl;
            exit(255);
        }
        if (options.count("fpim") && options.count("fpjm")) {
            std::cerr << "error parsing options: --fpim and --fpjm exclude each other" << std::endl;
            exit(255);
        }

        return options;
    } catch (const cxxopts::exceptions::exception & e) {
        // A bad argv must exit NON-zero (10/20 are verdicts, 255 is the error
        // convention): benchmark harnesses map exit codes, and exiting 0 here
        // silently turned option mistakes into fake "success" runs.
        std::cerr << "error parsing options: " << e.what() << std::endl;
        exit(255);
    }
}


VOID_TASK_0(gc_start)
{
    std::cerr << "starting garbage collection..." << std::endl;
}


VOID_TASK_0(gc_end)
{
    std::cerr << "garbage collection finished." << std::endl;
}


TASK_1(int, main_task, cxxopts::ParseResult*, _options)
{
    auto & options = *_options;

    bool verbose = options["verbose"].count() > 0;

    // First initialize the HOA data structure
    const auto t_before_parsing = wctime();
    auto data = (HoaData*)malloc(sizeof(HoaData));
    defaultsHoa(data);

    if (options.unmatched().empty()) {
        int ret = parseHoa(stdin, data);
        if (ret != 0) return ret;
    } else {
        std::string filename = options.unmatched()[0];
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
    if (verbose) {
        std::cerr << "\033[1;37mfinished parsing automaton in " << std::fixed << (t_after_parsing - t_before_parsing) << " sec.\033[m" << std::endl;
        std::cerr << "automaton has " << data->noStates << " states." << std::endl;
    }

    // First check if the automaton is a parity automaton (after rewriting
    // acceptance conditions that are parity conditions in disguise)
    bool isMaxParity = true;
    short controllerParity = 0;
    const int originalStates = data->noStates;
    int ret = paritizeHoa(data, options["paritize"].count() > 0);
    if (ret == 104) std::cerr << "(use --paritize to add the memory)" << std::endl;
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

    // Initialize Sylvan
    sylvan_set_limits(512LL << 22, 1, 14); // should be enough (2 gigabytes)
    sylvan_init_package();
    sylvan_init_mtbdd();
    sylvan_init_zdd();
    if (verbose) sylvan_gc_hook_pregc(TASK(gc_start));
    if (verbose) sylvan_gc_hook_postgc(TASK(gc_end));

    bool explicit_solver = options["sym"].count() == 0;
    bool naive_splitting = options["naive"].count() > 0;
    bool explicit_splitting = options["explicit"].count() > 0;
    bool write_pg = options["print-game"].count() > 0;
    bool no_solve = options["no-solve"].count() > 0;
    bool bisim_game = options["bisim"].count() > 0 or options["bisim-game"].count() > 0;
    bool bisim_sol = options["bisim"].count() > 0 or options["bisim-sol"].count() > 0;

    std::unique_ptr<SymGame> sym = nullptr;
    bool realizable; // not known yet
    MTBDD anyMoves = mtbdd_false;   // with --mealy-any: every move into the winning region
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
            if (verbose) {
                std::cerr << "\033[1;37mfinished constructing game in " << std::fixed << (t_after_splitting - t_before_splitting) << " sec.\033[m" << std::endl;
            }
        } else if (explicit_splitting) {
            const double t_before_splitting = wctime();
            game.reset(GameConstructor::constructGame(data, isMaxParity, controllerIsOdd));
            const double t_after_splitting = wctime();
            if (verbose) {
                std::cerr << "\033[1;37mfinished constructing game in " << std::fixed << (t_after_splitting - t_before_splitting) << " sec.\033[m" << std::endl;
            }
        } else {
            const double t_1 = wctime();
            vstart = 0; // always set to 0 by constructSymGame
            sym = GameConstructor::constructSymGame(data, isMaxParity, controllerIsOdd);
            const double t_2 = wctime();
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
                if (verbose) std::cerr << "\033[1;37mfinished bisimulation minimisation of game in " << std::fixed << (t_after - t_before) << " sec.\033[m" << std::endl;
            }
            const double t_3 = wctime();
            game = sym->toExplicit(vertex_to_bdd);
            const double t_4 = wctime();
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

        std::string solver = "tl";
        for (const auto & id : pg::Solvers::getSolverIDs()) {
            if (options.count(id)) solver = id;
        }

        pg::Oink engine(*game, verbose ? std::cerr : log);
        engine.setTrace(0); //verbose ? 1 : 0); actually don´t -- maybe add a 2nd verbosity level later
        engine.setRenumber();
        engine.setSolver(solver);
        engine.setWorkers(-1);

        // and run the solver
        engine.run();
        double end = wctime();

        // report how long it all took
        if (verbose) {
            std::cerr << "\033[1;37mfinished solving game in " << std::fixed << (end - begin) << " sec.\033[m" << std::endl;
        }

        // undo the sorting ## FOR SOME REASON, THIS CAN BE SLOW!
        game->permute(mapping.data());
        mapping.clear(); // erase result, we don't need it anymore

        realizable = game->getWinner(vstart) == 0;

        // finally, check if the initial vertex is won by controller or environment, and unless
        // only the verdict is asked for, get the strategy
        if (realizable && options["real"].count() == 0) {
            if (sym != nullptr && options["mealy-any"].count() > 0) {
                std::map<MTBDD, std::vector<MTBDD>> str;  // cap to winning priostates
                for (int v=0; v<game->vertexcount(); v++) {
                    if (game->owner(v) != 0 || game->getWinner(v) != 0) continue;
                    for (const int* e = game->outs(v); *e != -1; e++) {
                        if (game->getWinner(*e) == 0) str[vertex_to_bdd.at(v)].push_back(vertex_to_bdd.at(*e));
                    }
                }
                if (sym->applyMultiStrategy(str)) anyMoves = sym->strategies;
            }
            const bool fpim = options["fpim"].count() > 0;
            const bool fpjm = options["fpjm"].count() > 0;
            if (sym != nullptr && (fpim || fpjm || options["permissive"].count() > 0)) {
                auto allowed = PermissiveStrategy::positional(*game);
                if (fpim || fpjm) {
                    const double t_before = wctime();
                    const char* name = fpim ? "fixpoint iteration" : "fixpoint iteration with justifications";
                    std::vector<int> winner;
                    auto multi = fpim ? PermissiveStrategy::fpim(*game, winner) : PermissiveStrategy::fpjm(*game, winner);
                    bool agrees = true;
                    for (int v=0; v<game->vertexcount(); v++) agrees = agrees && winner[v] == game->getWinner(v);
                    if (agrees && PermissiveStrategy::verify(*game, vstart, multi)) {
                        allowed = std::move(multi);
                    } else {
                        // cannot happen; keep the positional strategy rather than risk a wrong controller
                        std::cerr << "error: " << name << " disagrees with " << solver << ", using its strategy" << std::endl;
                    }
                    if (verbose) {
                        long moves = 0;
                        for (const auto& a : allowed) moves += static_cast<long>(a.size());
                        std::cerr << name << ": " << moves << " moves of the controller" << std::endl;
                        std::cerr << "\033[1;37mfinished " << name << " in " << std::fixed << (wctime() - t_before) << " sec.\033[m" << std::endl;
                    }
                }
                std::vector<std::vector<int>> extended;  // the permissive extension of allowed, if any
                if (options["permissive"].count() > 0) {
                    const double t_before = wctime();
                    PermissiveStrategy::Stats stats;
                    extended = PermissiveStrategy::extend(*game, vstart, allowed, stats);
                    if (!PermissiveStrategy::verify(*game, vstart, extended)) {
                        // cannot happen; keep the strategy that was extended rather than risk a wrong controller
                        std::cerr << "error: the multi-strategy is not winning, using the strategy before extending it" << std::endl;
                        extended.clear();
                    }
                    if (verbose) {
                        std::cerr << "multi-strategy: " << stats.added << " of " << stats.candidates
                                  << " further edges allowed at " << stats.vertices << " reachable vertices" << std::endl;
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
                if (!extended.empty() && options["mealy-sat"].count() > 0) {
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
                    if (game->owner(v) == 0 && game->getWinner(v) == 0) {
                        auto a = vertex_to_bdd.at(v);
                        auto b = vertex_to_bdd.at(game->getStrategy(v));
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

        if (verbose) {
            std::cerr << "\033[1;37mfinished solving game in " << std::fixed << (t_after_solve - t_before_solve) << " sec.\033[m" << std::endl;
        }
    }

    if (options["real"].count() > 0) {
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

    if (options["mealy-sat"].count() > 0) {
        const auto t_before = wctime();
        MealyMinimizer::Options minimise;
        minimise.seconds = options["mealy-time"].as<double>();
        minimise.verbose = verbose;
        minimise.extended = extendedMoves;
        minimise.anyMoves = anyMoves;
        minimise.anySeconds = minimise.seconds;
        auto minimized = MealyMinimizer::minimize(*sym, minimise);
        if (minimized.game != nullptr) {
            sym = std::move(minimized.game);
            if (verbose) std::cerr << "minimised controller has " << minimized.states << " states." << std::endl;
        }
        if (verbose) {
            std::cerr << "\033[1;37mfinished mealy minimisation in " << std::fixed << (wctime() - t_before) << " sec.\033[m" << std::endl;
        }
    }
    mtbdd_unprotect(&anyMoves);
    mtbdd_unprotect(&extendedMoves);

    {
        const auto t_before = wctime();
        sym->postprocess(verbose);
        const auto t_after = wctime();

        if (verbose) {
            std::cerr << "\033[1;37mfinished post processing in " << std::fixed << (t_after - t_before) << " sec.\033[m" << std::endl;
        }
    }

    if (verbose) {
        // sym->print_trans();
        // sym->print_strategies();
    }

    const bool best = options["best"].count() > 0;

    if (bisim_sol || best) {
        const auto t_before = wctime();
        MTBDD partition = CALL(min_lts_strong, sym.get(), true);
        mtbdd_protect(&partition);
        if (verbose) {
            // CALL(print_partition, sym, partition);
        }
        CALL(minimize, sym.get(), partition, verbose);
        mtbdd_unprotect(&partition);
        const auto t_after = wctime();
        if (verbose) std::cerr << "\033[1;37mfinished bisimulation minimisation of solution in " << std::fixed << (t_after - t_before) << " sec.\033[m" << std::endl;
    }

    if (options["print-kiss"].count() > 0) {
        sym->print_kiss(true);
        if (verbose) {
            std::cerr << "\033[1;37mtotal time was " << std::fixed << (wctime() - t_before_parsing) << " sec.\033[m" << std::endl;
        }
        exit(10);
    }

    /**
     * maybe print witness parity game, which should be fully won by even
     */

    if (options["print-witness"].count() > 0) {
        auto pargame = sym->strategyToPG();
        pargame->write_pgsolver(std::cout);
        if (verbose) {
            std::cerr << "\033[1;37mtotal time was " << std::fixed << (wctime() - t_before_parsing) << " sec.\033[m" << std::endl;
        }
        exit(10);
    }

    if (best) {
        auto var1b = AIGEncoder(*data, *sym).setSop().encode();
        auto var2b = AIGEncoder(*data, *sym).setOneHot().setIsop().encode();
        auto var3b = AIGEncoder(*data, *sym).setOneHot().encode();
        auto var4b = AIGEncoder(*data, *sym).encode();

        if (verbose) {
            std::cerr << "bisim, sop: " << var1b->getNumAnds() << " gates and " << var1b->getNumLatches() << " latches" << std::endl;
            std::cerr << "bisim, onehot+isop: " << var2b->getNumAnds() << " gates and " << var2b->getNumLatches() << " latches" << std::endl;
            std::cerr << "bisim, onehot+ite: " << var3b->getNumAnds() << " gates and " << var3b->getNumLatches() << " latches" << std::endl;
            std::cerr << "bisim, binary+ite: " << var4b->getNumAnds() << " gates and " << var4b->getNumLatches() << " latches" << std::endl;
        }

        if (options["compress"].count() > 0) {
            ABCMinimization(*var1b, verbose).compress(compress_timeout/4);
            ABCMinimization(*var2b, verbose).compress(compress_timeout/4);
            ABCMinimization(*var3b, verbose).compress(compress_timeout/4);
            ABCMinimization(*var4b, verbose).compress(compress_timeout/4);

            if (verbose) {
                std::cerr << "sizes after compressing with ABC:" << std::endl;
                std::cerr << "bisim, sop: " << var1b->getNumAnds() << " gates and " << var1b->getNumLatches() << " latches" << std::endl;
                std::cerr << "bisim, onehot+isop: " << var2b->getNumAnds() << " gates and " << var2b->getNumLatches() << " latches" << std::endl;
                std::cerr << "bisim, onehot+ite: " << var3b->getNumAnds() << " gates and " << var3b->getNumLatches() << " latches" << std::endl;
                std::cerr << "bisim, binary+ite: " << var4b->getNumAnds() << " gates and " << var4b->getNumLatches() << " latches" << std::endl;
            }
        }

        auto var1bsize = var1b->getNumAnds() + var1b->getNumLatches();
        auto var2bsize = var2b->getNumAnds() + var2b->getNumLatches();
        auto var3bsize = var3b->getNumAnds() + var3b->getNumLatches();
        auto var4bsize = var4b->getNumAnds() + var4b->getNumLatches();

        auto smallest_size = var1bsize;
        auto * smallest = &var1b;
        if (smallest_size > var2bsize) {
            smallest_size = var2bsize;
            smallest = &var2b;
        }
        if (smallest_size > var3bsize) {
            smallest_size = var3bsize;
            smallest = &var3b;
        }
        if (smallest_size > var4bsize) {
            smallest_size = var4bsize;
            smallest = &var4b;
        }

        if (options.count("write-binary")) {
            (*smallest)->writeBinary(stdout);
        } else if (options["write-ascii"].count() > 0) {
            (*smallest)->writeAscii(stdout);
        }

        if (verbose) {
            std::cerr << "\033[1;37mtotal time was " << std::fixed << (wctime() - t_before_parsing) << " sec.\033[m" << std::endl;
        }
        exit(10);
    }

    AIGEncoder encoder(*data, *sym);
    if (verbose) encoder.setVerbose();
    if (options["isop"].count() > 0) encoder.setIsop();
    if (options["onehot"].count() > 0) encoder.setOneHot();
    if (options["sop"].count() > 0) encoder.setSop();
    const double t_before_enc = wctime();
    auto circuit = encoder.encode();
    const double t_after_enc = wctime();
    if (verbose) std::cerr << "\033[1;37mfinished encoding in " << std::fixed << (t_after_enc - t_before_enc) << " sec.\033[m" << std::endl;

    /**
     * maybe compress with ABC
     */
    if (options["drewrite"].count() > 0) {
        if (verbose) std::cerr << "size of AIG before drw+drf: " << circuit->getNumAnds() << " gates." << std::endl;
        const double t_before = wctime();
        ABCMinimization(*circuit, verbose).drewrite(compress_timeout);
        const double t_after = wctime();
        if (verbose) std::cerr << "size of AIG after drw+drf: " << circuit->getNumAnds() << " gates." << std::endl;
        if (verbose) std::cerr << "\033[1;37mfinished drw+drf in " << std::fixed << (t_after - t_before) << " sec.\033[m" << std::endl;
    }

    if (options["compress"].count() > 0) {
        if (verbose) std::cerr << "size of AIG before compression: " << circuit->getNumAnds() << " gates." << std::endl;
        const double t_before = wctime();
        ABCMinimization(*circuit, verbose).compress(compress_timeout);
        const double t_after = wctime();
        if (verbose) std::cerr << "size of AIG after compression: " << circuit->getNumAnds() << " gates." << std::endl;
        if (verbose) std::cerr << "\033[1;37mfinished compression in " << std::fixed << (t_after - t_before) << " sec.\033[m" << std::endl;
    }

    if (verbose) {
        std::cerr << "final size of AIG: " << circuit->getNumAnds() << " gates." << std::endl;
    }

    if (options["write-binary"].count() > 0) {
        circuit->writeBinary(stdout);
    } else if (options["write-ascii"].count() > 0) {
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
    auto options = handleOptions(argc, argv);

    // Initialize Lace, only 1 worker
    lace_start(1, 1024*1024*2); // initialize Lace, but sequentially
                                // also get a large enough task size... (2M tasks) for PSI!

    int res = RUN(main_task, &options);

    lace_stop();

    return res;
}
