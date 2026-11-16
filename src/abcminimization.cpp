/**
 * Copyright Tom van Dijk
 */

#include <abcminimization.hpp>
#include <algorithm>
#include <chrono>
#include <stdexcept>
#include <iostream>
#include <cstdio>
#include <filesystem>
#include <sstream>
#include <unistd.h>
#include <vector>

extern "C" {
#include "aig/gia/gia.h"
}

// slightly modified version of 'alias compress2rs' from 'abc.rc' file
const std::vector<std::string> ABCMinimization::compressCommands ({
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

namespace {

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

}  // namespace

void ABCMinimization::compress(int seconds, bool rrr)
{
    // without outputs no logic is observable; ABC has nothing to minimise and warns at every command
    if (circuit.getNumOutputs() == 0 || seconds <= 0) return;

    const auto end = std::chrono::steady_clock::now() + std::chrono::seconds(seconds);
    auto left = [&]() -> long {   // milliseconds
        return std::chrono::duration_cast<std::chrono::milliseconds>(end - std::chrono::steady_clock::now()).count();
    };

    AbcSession session;
    Abc_Frame_t* pAbc = Abc_FrameGetGlobalFrame();

    writeToAbc(pAbc);
    const int latches = Abc_NtkLatchNum(Abc_FrameReadNtk(pAbc));

    // Rounds of the compress commands and, with rrr, resubstitution with redundancy analysis (&rrr;
    // Miyasaka, Mishchenko, Wawrzynek, IWLS 2025), multi-add (-X 1: deterministic, and much faster than
    // single-add on the large circuits). &rrr works on the combinational logic and drops the registers (the
    // latch outputs become inputs, the latch inputs outputs, both last), which &setregnum restores. It runs on
    // a copy, within the time that is left (when the time is up it keeps what it found), and the copy replaces
    // the network only if it is smaller: a round that runs out of time, or an &rrr that fails, costs only its
    // time. The rounds go on while they gain something, with &rrr while they gain at least 2%: its rounds cost
    // minutes on the largest circuits and the later ones gain little. Stops when the time is up (also between
    // the compress commands), and before a round that would take longer than the time that is left.
    int size = getAbcNetworkSize(pAbc);
    while (size > 0) {
        const long before = left();
        executeCompressCommands(pAbc, end);
        if (rrr && left() > 0) {
            try {
                executeAbcCommand(pAbc, "&get -n");
                executeAbcCommand(pAbc, "&rrr -X 1 -T " + std::to_string(std::max(1L, left() / 1000)));
                executeAbcCommand(pAbc, "&setregnum -N " + std::to_string(latches));
                Gia_Man_t* copy = Abc_FrameReadGia(pAbc);
                if (copy != nullptr && Gia_ManAndNum(copy) < getAbcNetworkSize(pAbc)) executeAbcCommand(pAbc, "&put");
            } catch (const std::runtime_error&) {
                // &rrr does not take every network (e.g. one without inputs): the network stays as it is
            }
        }
        const int next = getAbcNetworkSize(pAbc);
        const bool gained = rrr ? 50L * (size - next) >= size : next < size;
        size = next;
        const long after = left();
        if (!gained || before - after >= after) break;
    }

    readFromAbc(pAbc);
}

void ABCMinimization::drewrite(int seconds)
{
    if (circuit.getNumOutputs() == 0 || seconds <= 0) return;

    auto start = std::chrono::steady_clock::now();
    AbcSession session;
    Abc_Frame_t* pAbc = Abc_FrameGetGlobalFrame();

    writeToAbc(pAbc);

    // compress until convergence
    int new_num_nodes = getAbcNetworkSize(pAbc);
    int old_num_nodes = new_num_nodes + 1;
    while (new_num_nodes > 0 && new_num_nodes < old_num_nodes) {
        auto start2 = std::chrono::steady_clock::now();
        executeAbcCommand(pAbc, "drw");
        executeAbcCommand(pAbc, "balance");
        executeAbcCommand(pAbc, "drf");
        executeAbcCommand(pAbc, "dc2");
        old_num_nodes = new_num_nodes;
        new_num_nodes = getAbcNetworkSize(pAbc);
        if ((old_num_nodes-new_num_nodes)<old_num_nodes/100) break; // 1% improvement or better pls
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start2).count();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
        if ((elapsed+duration) >= 1000L * seconds) break;
    }

    readFromAbc(pAbc);
}

void ABCMinimization::executeAbcCommand(Abc_Frame_t* pAbc, const std::string& command) const {
    if (Cmd_CommandExecute( pAbc, command.c_str())) {
        throw std::runtime_error("Cannot execute ABC command: " + command);
    }
    if (verbose) std::cerr << "after " << command << ": " << getAbcNetworkSize(pAbc) << std::endl;
}

void ABCMinimization::executeCompressCommands(Abc_Frame_t* pAbc, std::chrono::steady_clock::time_point end) const {
    for (const auto& command : compressCommands) {
        if (std::chrono::steady_clock::now() >= end) return;
        executeAbcCommand(pAbc, command);
    }
}

int ABCMinimization::getAbcNetworkSize(Abc_Frame_t* pAbc) {
    Abc_Ntk_t* pNtk = Abc_FrameReadNtk(pAbc);
    return Abc_NtkNodeNum(pNtk);
}

// A new temporary file in the system's temporary directory (not the working directory, which need
// not be writable): its descriptor, and its name in filename.
int ABCMinimization::getTmpFile(std::string& filename) {
    std::string pattern = (std::filesystem::temp_directory_path() / "knor.XXXXXX").string();
    std::vector<char> name(pattern.begin(), pattern.end());
    name.push_back('\0');
    int fd = mkstemp(name.data());
    if (fd == -1) {
        throw std::runtime_error("Could not create temporary file: " + pattern);
    }
    filename = name.data();
    return fd;
}

void ABCMinimization::writeToAbc(Abc_Frame_t* pAbc) const {
    std::string tmp_filename;
    int fd = getTmpFile(tmp_filename);

    // write AIGER out to be read by ABC
    FILE* file = fdopen(fd, "w");
    if (file == nullptr) {
        throw std::runtime_error("Could not open temporary file: " + tmp_filename);
    }
    int write_result = circuit.writeBinary(file);
    fclose(file);
    if (write_result == 0) {
        throw std::runtime_error("Could not write AIGER circuit to file: " + tmp_filename);
    }

    std::stringstream read_command;
    read_command << "read_aiger " << tmp_filename;
    executeAbcCommand(pAbc, read_command.str());

    std::remove(tmp_filename.c_str());
}

void ABCMinimization::readFromAbc(Abc_Frame_t* pAbc) {
    std::string tmp_filename;
    int fd = getTmpFile(tmp_filename);

    std::stringstream write_command;
    write_command << "write_aiger -s " << tmp_filename;
    executeAbcCommand(pAbc, write_command.str());

    // read AIGER back, delete comments added by ABC
    FILE* file = fdopen(fd, "r");
    if (file == nullptr) {
        throw std::runtime_error("Could not open temporary file: " + tmp_filename);
    }
    // read_aiger
    circuit.readFile(file);
    fclose(file);
    std::remove(tmp_filename.c_str());
}
