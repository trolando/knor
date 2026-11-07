/**
 * Copyright Tom van Dijk
 */

#include <abcminimization.hpp>
#include <chrono>
#include <stdexcept>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <sstream>
#include <unistd.h>

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

void ABCMinimization::compress(int timeout)
{
    // without outputs no logic is observable; ABC has nothing to minimise and warns at every command
    if (circuit.getNumOutputs() == 0) return;

    auto start = std::chrono::steady_clock::now();

    AbcSession session;
    Abc_Frame_t* pAbc = Abc_FrameGetGlobalFrame();

    writeToAbc(pAbc);

    // compress until convergence
    int new_num_nodes = getAbcNetworkSize(pAbc);
    int old_num_nodes = new_num_nodes + 1;
    while (new_num_nodes > 0 && new_num_nodes < old_num_nodes) {
        auto start2 = std::chrono::steady_clock::now();
        executeCompressCommands(pAbc);
        old_num_nodes = new_num_nodes;
        new_num_nodes = getAbcNetworkSize(pAbc);
        // if ((old_num_nodes-new_num_nodes)<old_num_nodes/100) break; // 2.5% improvement or better pls
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start2).count();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
        if ((elapsed+duration) >= timeout) break;
    }

    readFromAbc(pAbc);
}

void ABCMinimization::drewrite(int timeout)
{
    if (circuit.getNumOutputs() == 0) return;

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
        if ((elapsed+duration) >= timeout) break;
    }

    readFromAbc(pAbc);
}

void ABCMinimization::executeAbcCommand(Abc_Frame_t* pAbc, const std::string& command) const {
    if (Cmd_CommandExecute( pAbc, command.c_str())) {
        throw std::runtime_error("Cannot execute ABC command: " + command);
    }
    if (verbose) std::cerr << "after " << command << ": " << getAbcNetworkSize(pAbc) << std::endl;
}

void ABCMinimization::executeCompressCommands(Abc_Frame_t* pAbc) const {
    for (const auto& command : compressCommands) {
        executeAbcCommand(pAbc, command);
    }
}

int ABCMinimization::getAbcNetworkSize(Abc_Frame_t* pAbc) {
    Abc_Ntk_t* pNtk = Abc_FrameReadNtk(pAbc);
    return Abc_NtkNodeNum(pNtk);
}

int ABCMinimization::getTmpFile(char* tmp_filename) {
    // TODO this needs to be fixed (maybe to std::string?) strcpy usage is technically prone to errors
    std::strcpy(tmp_filename, "knor.XXXXXX");
    int fd = mkstemp(tmp_filename);
    if (fd == -1) {
        throw std::runtime_error("Could not create temporary file: " + std::string(tmp_filename));
    }
    return fd;
}

void ABCMinimization::writeToAbc(Abc_Frame_t* pAbc) const {
    char tmp_filename[256];
    int fd = getTmpFile(tmp_filename);

    // write AIGER out to be read by ABC
    FILE* file = fdopen(fd, "w");
    if (file == nullptr) {
        throw std::runtime_error("Could not open temporary file: " + std::string(tmp_filename));
    }
    int write_result = circuit.writeBinary(file);
    fclose(file);
    if (write_result == 0) {
        throw std::runtime_error("Could not write AIGER circuit to file: " + std::string(tmp_filename));
    }

    std::stringstream read_command;
    read_command << "read_aiger " << tmp_filename;
    executeAbcCommand(pAbc, read_command.str());

    std::remove(tmp_filename);
}

void ABCMinimization::readFromAbc(Abc_Frame_t* pAbc) {
    char tmp_filename[256];
    int fd = getTmpFile(tmp_filename);

    std::stringstream write_command;
    write_command << "write_aiger -s " << tmp_filename;
    executeAbcCommand(pAbc, write_command.str());

    // read AIGER back, delete comments added by ABC
    FILE* file = fdopen(fd, "r");
    if (file == nullptr) {
        throw std::runtime_error("Could not open temporary file: " + std::string(tmp_filename));
    }
    // read_aiger
    circuit.readFile(file);
    fclose(file);
    std::remove(tmp_filename);
}
