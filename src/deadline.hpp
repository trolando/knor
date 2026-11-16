/**
 * Copyright Tom van Dijk
 */

#ifndef KNOR_DEADLINE_HPP
#define KNOR_DEADLINE_HPP

#include <algorithm>
#include <chrono>

/**
 * The end of the time of the run (--timeout). The stages that can stop early with what they have found
 * stop when it has passed: the permissive extension, the SAT searches of the minimiser, the encoding
 * search and the compression. Constructing and solving the game cannot; they are not limited.
 */
class Deadline final {
public:
    Deadline() = default;   // none
    explicit Deadline(int seconds)
        : limited(true), end(std::chrono::steady_clock::now() + std::chrono::seconds(seconds)) {}

    bool passed() const { return limited && std::chrono::steady_clock::now() >= end; }

    // The whole seconds that are left, at most limit (without a deadline, limit).
    int seconds(int limit) const
    {
        if (!limited) return limit;
        const auto left = std::chrono::duration_cast<std::chrono::seconds>(end - std::chrono::steady_clock::now());
        return static_cast<int>(std::clamp<long long>(left.count(), 0, limit));
    }

private:
    bool limited = false;
    std::chrono::steady_clock::time_point end;
};

#endif
