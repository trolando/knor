/**
 * Copyright Tom van Dijk
 */

// Tests StateSearch, the search over the number of controller states of the
// minimiser, against an oracle instead of the SAT encoding: the specification
// needs `minimum` states, a call with at least that many finds a controller
// with between `minimum` and k reachable states, a call with fewer refutes k,
// except for the numbers that the oracle leaves undecided and the encodings
// that are too large (from some number on). Every run must terminate, try each
// number at most once and only between the lower bound and the known
// controller, and conclude only what the oracle allows. With fewerFirst, it
// must go below an undecided number when it can and stop only when nothing is
// left to try (without, the search goes upwards past the undecided numbers and
// leaves the numbers it stepped over). Each search also runs with unit steps,
// which stop at the first undecided number.
// Exhaustive for small searches, random (with a fixed seed) for larger ones.

#include <cstdint>
#include <iostream>
#include <set>
#include <sstream>
#include <string>

#include <statesearch.hpp>

namespace {

struct Scenario {
    int lowerBound;
    int known;              // the states of the known controller (at least minimum)
    int minimum;            // the true minimum, at least lowerBound (which is at most known)
    std::set<int> undecided;
    int tooLargeFrom;       // the encodings for this many states or more are too large (0: none)
    uint32_t seed;          // chooses the reachable states of the controllers found
    bool fewerFirst;
    bool unit = false;
};

std::string describe(const Scenario& sc)
{
    std::ostringstream out;
    out << "lowerBound " << sc.lowerBound << ", known " << sc.known << ", minimum " << sc.minimum
        << ", undecided {";
    for (int u : sc.undecided) out << " " << u;
    out << " }, too large from " << sc.tooLargeFrom << ", seed " << sc.seed
        << (sc.fewerFirst ? ", fewer first" : "") << (sc.unit ? ", unit steps" : "");
    return out.str();
}

uint32_t nextRandom(uint32_t& state)
{
    state = state * 1664525u + 1013904223u;
    return state >> 8;
}

int failures = 0;

void fail(const Scenario& sc, const std::string& what)
{
    if (failures++ < 20) std::cerr << "FAIL: " << what << " (" << describe(sc) << ")" << std::endl;
}

// Runs the search as the minimiser does and checks it; returns the number of calls.
int run(const Scenario& sc)
{
    StateSearch sizes(sc.lowerBound, sc.known, sc.fewerFirst, sc.unit);
    uint32_t random = sc.seed;
    std::set<int> tried, undecidedSeen;
    int refuted = std::max(1, sc.lowerBound) - 1;    // what the oracle has refuted so far
    int found = sc.known;
    int tooLarge = sc.known;                          // the smallest number found too large
    int calls = 0;
    for (int k = sizes.next(); k > 0; k = sizes.next()) {
        if (++calls > 4 * sc.known + 16) {
            fail(sc, "does not terminate");
            return calls;
        }
        if (k < sc.lowerBound || k >= sc.known) fail(sc, "tries " + std::to_string(k) + " states, out of range");
        if (!tried.insert(k).second) fail(sc, "tries " + std::to_string(k) + " states twice");
        if (k >= found || k >= tooLarge) fail(sc, "tries " + std::to_string(k) + " states, not below the best");
        if (k <= refuted) fail(sc, "tries " + std::to_string(k) + " states, already refuted");
        if (sc.tooLargeFrom > 0 && k >= sc.tooLargeFrom) {
            sizes.tooLarge(k);
            tooLarge = std::min(tooLarge, k);
        } else if (sc.undecided.count(k)) {
            const bool open = [&] {
                for (int j = refuted + 1; j < k; j++) {
                    if (!undecidedSeen.count(j)) return true;
                }
                return false;
            }();
            sizes.undecide(k);
            undecidedSeen.insert(k);
            // fewerFirst: a number below k that is still open comes next
            const int after = sizes.next();
            if (sc.fewerFirst && open && (after == 0 || after > k)) {
                fail(sc, "after undecided " + std::to_string(k) + " tries " + std::to_string(after) + ", not fewer");
            }
        } else if (k >= sc.minimum) {
            const int states = sc.minimum + static_cast<int>(nextRandom(random) % (k - sc.minimum + 1));
            sizes.suffice(states);
            found = std::min(found, states);
        } else {
            sizes.refute(k);
            refuted = std::max(refuted, k);
        }
    }

    if (sizes.states() != found) fail(sc, "states " + std::to_string(sizes.states()) + ", expected " +
                                          std::to_string(found));
    if (sizes.proven() && sc.minimum < found) fail(sc, "proven, but " + std::to_string(found) + " is not minimal");
    // with fewerFirst, it stops only when every number between the refuted ones and the best is
    // undecided (or too large); unit steps also stop above the first undecided number, and only
    // there: every number below it is refuted
    for (int j = refuted + 1; sc.fewerFirst && j < std::min(found, tooLarge); j++) {
        if (sc.unit && !undecidedSeen.empty() && j > *undecidedSeen.begin()) continue;
        if (!undecidedSeen.count(j)) fail(sc, "stops with " + std::to_string(j) + " states untried");
    }
    if (sc.unit && !undecidedSeen.empty() && *undecidedSeen.begin() != refuted + 1 && tooLarge == sc.known
        && found == sc.known) {
        fail(sc, "unit steps leave a number below the first undecided one");
    }
    // with an oracle that decides everything, the search finds and proves the minimum
    if (sc.undecided.empty() && sc.tooLargeFrom == 0 && sc.lowerBound <= sc.minimum
        && (found != sc.minimum || !sizes.proven())) {
        fail(sc, "does not find and prove the minimum");
    }
    return calls;
}

}  // namespace

int main()
{
    long scenarios = 0, calls = 0;
    // exhaustive: every lower bound, known controller, minimum and set of undecided numbers
    for (int known = 1; known <= 9; known++) {
        for (int lowerBound = 0; lowerBound <= known; lowerBound++) {
            for (int minimum = std::max(1, lowerBound); minimum <= known; minimum++) {
                const int first = std::max(1, lowerBound);
                const int range = std::max(0, known - first);
                for (uint32_t mask = 0; mask < (1u << range); mask++) {
                    std::set<int> undecided;
                    for (int i = 0; i < range; i++) {
                        if (mask & (1u << i)) undecided.insert(first + i);
                    }
                    for (int tooLargeFrom : {0, first + range / 2 + 1}) {
                        for (bool fewerFirst : {false, true}) {
                            for (uint32_t seed : {1u, 2u, 3u}) {
                                for (bool unit : {false, true}) {
                                    calls += run({lowerBound, known, minimum, undecided, tooLargeFrom, seed, fewerFirst,
                                                  unit});
                                    scenarios++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // random: larger searches
    uint32_t random = 12345;
    for (int i = 0; i < 20000; i++) {
        const int known = 2 + static_cast<int>(nextRandom(random) % 200);
        const int lowerBound = static_cast<int>(nextRandom(random) % (known + 1));
        const int first = std::max(1, lowerBound);
        const int minimum = first + static_cast<int>(nextRandom(random) % (known - first + 1));
        std::set<int> undecided;
        const int density = static_cast<int>(nextRandom(random) % 4);
        for (int k = first; k < known; k++) {
            if (density > 0 && nextRandom(random) % (2 * density + 2) == 0) undecided.insert(k);
        }
        const int tooLargeFrom = nextRandom(random) % 3 == 0 ? first + static_cast<int>(nextRandom(random) % known) : 0;
        const uint32_t seed = nextRandom(random);
        const bool fewerFirst = nextRandom(random) % 2 == 0;
        calls += run({lowerBound, known, minimum, undecided, tooLargeFrom, seed, fewerFirst, nextRandom(random) % 2 == 0});
        scenarios++;
    }
    if (failures > 0) {
        std::cerr << failures << " failures in " << scenarios << " searches" << std::endl;
        return 1;
    }
    std::cout << scenarios << " searches (" << calls << " calls) as expected" << std::endl;
    return 0;
}
