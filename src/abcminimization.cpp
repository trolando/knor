/**
 * Copyright Tom van Dijk
 */

#include <abcminimization.hpp>
#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <functional>
#include <iostream>
#include <signal.h>
#include <spawn.h>
#include <stdexcept>
#include <sys/wait.h>
#include <thread>
#include <unistd.h>
#include <vector>
#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

extern "C" {
#include "misc/util/abc_namespaces.h"
#include "base/main/abcapis.h"
#include "base/abc/abc.h"
#include "base/main/main.h"
#include "aig/gia/gia.h"
}

extern char** environ;

/*
 * ABC checks a time limit only between the steps of a command, and a step can take far longer than the time that
 * is left: &rrr's redundancy removal checks no time at all, and with a limit of a minute ran for more than 15
 * minutes on a GR(1) controller of 15,000 gates. So the rounds run in a worker: knor itself, started anew (not
 * forked: a forked process shares the memory of knor, and counts it again in the memory of the process tree). It
 * reads the circuit from a file, writes each smaller network it finds to another (a new file renamed into place:
 * a stopped worker leaves a whole one), and is stopped when the time is up; the circuit becomes the last network
 * written, also when a round fails (an error, or ABC crashing). Without a worker (the executable is not known, or
 * cannot be started) the rounds run in this process, stopped only by their own checks of the time.
 */

namespace {

using Clock = std::chrono::steady_clock;

// slightly modified version of 'alias compress2rs' from 'abc.rc' file
const std::vector<std::string> COMPRESS_COMMANDS({
       "dc2",
       "balance -l",
       "resub -K 6 -l",
       "rewrite -l",
       "resub -K 6 -N 2",
       "refactor -l",
       "resub -K 8 -l",
       "balance -l",
       "resub -K 8 -N 2 -l",
       "rewrite -l",
       "resub -K 10 -l",
       "rewrite -z -l",
       "resub -K 10 -N 2 -l",
       "balance -l",
       "resub -K 12 -l",
       "refactor -z -l",
       "resub -K 12 -N 2 -l",
       "balance -l",
       "rewrite -z -l",
       "balance -l",
       "drw -z -C 100 -N 10000",
       "drf -z -C 10000 -K 15",
       "ifraig -C 20",
});

/**
 * Runs ABC for the lifetime of the object. ABC prints its messages to stdout, where knor writes the
 * circuit, so stdout goes to stderr meanwhile. ABC also keeps a command history in abc.history in
 * the working directory: batch mode records no commands, and the file is removed again unless it
 * was there before.
 */
class AbcSession final {
public:
    AbcSession() : had_history(std::filesystem::exists(HISTORY, error))
    {
        std::fflush(stdout);
        saved_stdout = dup(STDOUT_FILENO);
        if (saved_stdout != -1) dup2(STDERR_FILENO, STDOUT_FILENO);
        Abc_Start();
        Abc_FrameSetBatchMode(1);
    }

    ~AbcSession()
    {
        Abc_Stop();
        std::fflush(stdout);
        if (saved_stdout != -1) {
            dup2(saved_stdout, STDOUT_FILENO);
            close(saved_stdout);
        }
        if (!had_history) std::filesystem::remove(HISTORY, error);
    }

    AbcSession(const AbcSession&) = delete;
    AbcSession& operator=(const AbcSession&) = delete;

private:
    static constexpr const char* HISTORY = "abc.history";
    std::error_code error;
    bool had_history;
    int saved_stdout;
};

/**
 * The rounds of ABC commands on the circuit in the input file. After every command, a network with fewer gates
 * than the circuit so far is written to the result file.
 */
class Rounds final {
public:
    Rounds(const std::string& input, const std::string& result, long kept, bool verbose)
        : result(result), part(result + ".part"), kept(kept), verbose(verbose)
    {
        command("read_aiger " + input);
    }

    // Rounds of the compress commands while they gain: they are cheap, and one &rrr can take all the time that is
    // left. Then, with rrr, rounds of resubstitution with redundancy analysis (&rrr; Miyasaka, Mishchenko,
    // Wawrzynek, IWLS 2025) and the compress commands while they gain at least 2%: the rounds of &rrr cost minutes
    // on the largest circuits, and the later ones gain little. &rrr runs multi-add (-X 1: deterministic, and much
    // faster than single-add on the large circuits) on the combinational logic, and drops the registers (the latch
    // outputs become inputs, the latch inputs outputs, both last), which &setregnum restores. It runs on a copy,
    // within 90% of the time that is left (it checks the time only between its steps, and it keeps what it found
    // only if it stops itself), and the copy replaces the network only if it is smaller: an &rrr that runs out of
    // time, or that fails, costs only its time. Stops when the time is up (also between the compress commands).
    void compress(Clock::time_point end, bool rrr)
    {
        rounds(end, [&]() { compressCommands(end); }, [](int current, int next) { return next < current; });
        if (!rrr) return;
        const int latches = Abc_NtkLatchNum(Abc_FrameReadNtk(pAbc));
        rounds(end, [&]() {
            resubstitute(end, latches);
            compressCommands(end);
        }, [](int current, int next) { return 50L * (current - next) >= current; });
    }

    // rounds of drw, balance, drf and dc2 while they gain at least 1%
    void drewrite(Clock::time_point end)
    {
        rounds(end, [&]() {
            for (const char* c : {"drw", "balance", "drf", "dc2"}) {
                command(c);
                keep();
            }
        }, [](int current, int next) { return next < current && current - next >= current / 100; });
    }

private:
    AbcSession session;
    Abc_Frame_t* pAbc = Abc_FrameGetGlobalFrame();
    const std::string result;
    const std::string part;
    long kept;
    const bool verbose;

    // the round while it gains, and while it fits in the time that is left (it takes as long as the last time)
    void rounds(Clock::time_point end, const std::function<void()>& round, bool (*gained)(int current, int next))
    {
        int current = size();
        while (current > 0) {
            const auto start = Clock::now();
            round();
            const int next = size();
            const bool gain = gained(current, next);
            current = next;
            const auto now = Clock::now();
            if (!gain || now + (now - start) >= end) break;
        }
    }

    void compressCommands(Clock::time_point end)
    {
        for (const auto& c : COMPRESS_COMMANDS) {
            if (Clock::now() >= end) return;
            command(c);
            keep();
        }
    }

    void resubstitute(Clock::time_point end, int latches)
    {
        const long left = std::chrono::duration_cast<std::chrono::milliseconds>(end - Clock::now()).count();
        if (left <= 0) return;
        try {
            command("&get -n");
            command("&rrr -X 1 -T " + std::to_string(std::max(1L, left * 9 / 10000)));
            command("&setregnum -N " + std::to_string(latches));
            Gia_Man_t* copy = Abc_FrameReadGia(pAbc);
            if (copy != nullptr && Gia_ManAndNum(copy) < size()) {
                command("&put");
                keep();
            }
        } catch (const std::runtime_error&) {
            // &rrr does not take every network (e.g. one without inputs): the network stays as it is
        }
    }

    void command(const std::string& command) const
    {
        if (Cmd_CommandExecute(pAbc, command.c_str())) {
            throw std::runtime_error("Cannot execute ABC command: " + command);
        }
        if (verbose) std::cerr << "after " << command << ": " << size() << std::endl;
    }

    int size() const
    {
        return Abc_NtkNodeNum(Abc_FrameReadNtk(pAbc));
    }

    void keep()
    {
        const int gates = size();
        if (gates >= kept) return;
        std::remove(part.c_str());   // write_aiger reports no failure to write
        command("write_aiger -s " + part);
        if (std::rename(part.c_str(), result.c_str()) != 0) {
            throw std::runtime_error("Could not rename " + part + " to " + result);
        }
        kept = gates;
    }
};

// the rounds (drewrite, compress or compress-rrr) on the circuit in the input file, until the end
void work(const std::string& rounds, Clock::time_point end, const std::string& input, const std::string& result,
          long kept, bool verbose)
{
    Rounds r(input, result, kept, verbose);
    if (rounds == "drewrite") {
        r.drewrite(end);
    } else if (rounds == "compress" || rounds == "compress-rrr") {
        r.compress(end, rounds == "compress-rrr");
    } else {
        throw std::runtime_error("Unknown rounds: " + rounds);
    }
}

/**
 * A new directory in the system's temporary directory (not the working directory, which need not be writable),
 * removed with its files at the end of its scope.
 */
class TemporaryDirectory final {
public:
    TemporaryDirectory()
    {
        const std::string pattern = (std::filesystem::temp_directory_path() / "knor.XXXXXX").string();
        std::vector<char> name(pattern.begin(), pattern.end());
        name.push_back('\0');
        if (mkdtemp(name.data()) == nullptr) {
            throw std::runtime_error("Could not create temporary directory: " + pattern);
        }
        path = name.data();
    }

    ~TemporaryDirectory()
    {
        std::error_code error;
        std::filesystem::remove_all(path, error);
    }

    TemporaryDirectory(const TemporaryDirectory&) = delete;
    TemporaryDirectory& operator=(const TemporaryDirectory&) = delete;

    std::string file(const char* name) const
    {
        return (path / name).string();
    }

private:
    std::filesystem::path path;
};

// the path of the running executable, or empty if it is not known
std::string executable()
{
    std::error_code error;
#if defined(__linux__)
    const auto path = std::filesystem::read_symlink("/proc/self/exe", error);
#elif defined(__APPLE__)
    uint32_t length = 0;
    _NSGetExecutablePath(nullptr, &length);
    std::vector<char> name(length + 1, '\0');
    if (_NSGetExecutablePath(name.data(), &length) != 0) return "";
    const auto path = std::filesystem::canonical(name.data(), error);
#else
    const std::filesystem::path path;
#endif
    return error ? "" : path.string();
}

// Starts the executable with the arguments, with its stdout on stderr; its process, or -1 if it cannot be started.
pid_t spawn(const std::string& executable, const std::vector<std::string>& arguments)
{
    std::vector<char*> argv{const_cast<char*>(executable.c_str())};
    for (const auto& argument : arguments) argv.push_back(const_cast<char*>(argument.c_str()));
    argv.push_back(nullptr);
    posix_spawn_file_actions_t actions;
    if (posix_spawn_file_actions_init(&actions) != 0) return -1;
    pid_t pid = -1;
    const bool started = posix_spawn_file_actions_adddup2(&actions, STDERR_FILENO, STDOUT_FILENO) == 0
                         && posix_spawn(&pid, executable.c_str(), &actions, nullptr, argv.data(), environ) == 0;
    posix_spawn_file_actions_destroy(&actions);
    return started ? pid : -1;
}

// Waits for the process until the end, and stops it then.
void waitUntil(pid_t process, Clock::time_point end)
{
    for (;;) {
        const pid_t done = waitpid(process, nullptr, WNOHANG);
        if (done == process || (done == -1 && errno != EINTR)) return;
        if (Clock::now() >= end) {
            kill(process, SIGKILL);
            while (waitpid(process, nullptr, 0) == -1 && errno == EINTR) {}
            return;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
}

}  // namespace

void ABCMinimization::compress(int seconds, bool rrr)
{
    minimise(seconds, rrr ? "compress-rrr" : "compress");
}

void ABCMinimization::drewrite(int seconds)
{
    minimise(seconds, "drewrite");
}

void ABCMinimization::minimise(int seconds, const std::string& rounds)
{
    // without outputs no logic is observable; ABC has nothing to minimise and warns at every command
    if (circuit.getNumOutputs() == 0 || seconds <= 0) return;

    const auto end = Clock::now() + std::chrono::seconds(seconds);
    const TemporaryDirectory directory;
    const std::string input = directory.file("circuit.aig");
    const std::string result = directory.file("result.aig");   // the network written last, if any
    FILE* file = std::fopen(input.c_str(), "wb");
    if (file == nullptr) throw std::runtime_error("Could not open temporary file: " + input);
    const int written = circuit.writeBinary(file);
    std::fclose(file);
    if (written == 0) throw std::runtime_error("Could not write AIGER circuit to file: " + input);

    const long kept = circuit.getNumAnds();
    const std::string path = executable();
    const long milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - Clock::now()).count();
    const pid_t worker = path.empty() ? -1 : spawn(path, {WORKER, rounds, std::to_string(milliseconds), input,
                                                          result, std::to_string(kept), verbose ? "1" : "0"});
    if (worker == -1) {
        try {
            work(rounds, end, input, result, kept, verbose);
        } catch (const std::exception& e) {
            std::cerr << "compression stopped: " << e.what() << std::endl;
        }
    } else {
        waitUntil(worker, end);
    }

    file = std::fopen(result.c_str(), "rb");
    if (file == nullptr) return;   // nothing smaller found
    try {
        circuit.readFile(file);
    } catch (const std::runtime_error& e) {
        std::cerr << "compression lost: " << e.what() << std::endl;
    }
    std::fclose(file);
}

int ABCMinimization::worker(int argc, char* argv[])
{
    // knor WORKER ROUNDS MILLISECONDS INPUT RESULT GATES VERBOSE
    if (argc != 8) {
        std::cerr << "knor: " << WORKER << " takes the rounds, the milliseconds, the input, the result, the gates"
                  << " and the verbosity" << std::endl;
        return 255;
    }
    try {
        const long milliseconds = std::stol(argv[3]);
        // the worker stops itself as well, also when knor is gone
        signal(SIGALRM, SIG_DFL);
        alarm(static_cast<unsigned>(std::max(0L, milliseconds) / 1000 + 2));
        work(argv[2], Clock::now() + std::chrono::milliseconds(milliseconds), argv[4], argv[5], std::stol(argv[6]),
             std::strcmp(argv[7], "1") == 0);
    } catch (const std::exception& e) {
        std::cerr << "compression stopped: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
