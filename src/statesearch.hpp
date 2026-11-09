/**
 * Copyright Tom van Dijk
 */

#ifndef KNOR_STATESEARCH_HPP
#define KNOR_STATESEARCH_HPP

#include <algorithm>
#include <set>

// The search for the smallest number of states that suffice, from a lower bound
// to the states of a known controller: fewer states than a refuted number do not
// suffice either, and more than a number that suffices do. Upwards with growing
// steps (small encodings are much easier, so the search stays near the lower
// bound), and by bisection once a controller is found or an encoding is too large.
// With fewerFirst, after an undecided call first the numbers below it (in the
// middle of the lowest gap between the refuted and the undecided numbers): their
// encodings are smaller, and a controller with fewer states is often found much
// faster. With unit, upwards by one state, up to the first undecided number: where
// every further state multiplies the cost of a call (as observed when the
// candidates are checked for losing cycles), a call above the minimum costs more
// than all the calls below it, and a call above an undecided one, with half of
// its allowance, is unlikely to be decided (a scheduling choice: an undecided
// number never counts as refuted).
class StateSearch final {
public:
    StateSearch(int lowerBound, int known, bool fewerFirst, bool unit = false)
        : refuted(std::max(1, lowerBound) - 1), highest(refuted), found(known), top(known), fewerFirst(fewerFirst),
          unit(unit)
    {
    }

    // The next number of states to try, or 0 if every number below the bound is decided or undecided.
    int next() const
    {
        int low = refuted;
        for (int u : undecided) {
            if (fewerFirst && low + 1 < u) return low + (u - low) / 2;
            low = u;
        }
        if (low + 1 >= top) return 0;
        if (bisect) return low + (top - low) / 2;
        if (unit) return undecided.empty() ? low + 1 : 0;
        return std::min(low + step, top - 1);
    }

    // A controller with this many (reachable) states was found.
    void suffice(int states)
    {
        found = std::min(found, states);
        bound(found);
    }

    // The encoding for k states is too large: search below k.
    void tooLarge(int k) { bound(k); }

    void refute(int k)
    {
        refuted = std::max(refuted, k);
        undecided.erase(undecided.begin(), undecided.upper_bound(refuted));
        failedAt(k);
    }

    void undecide(int k)
    {
        undecided.insert(k);
        failedAt(k);
    }

    int states() const { return found; }
    bool proven() const { return refuted + 1 >= found; }

private:
    void bound(int k)
    {
        top = std::min(top, k);
        undecided.erase(undecided.lower_bound(top), undecided.end());
        bisect = true;
    }

    void failedAt(int k)
    {
        if (k > highest) {
            highest = k;
            if (!unit && step < MAX_STEP) step *= 2;
        }
    }

    static constexpr int MAX_STEP = 1 << 20;

    int refuted;                // no number up to it suffices
    int highest;                // the largest number that did not suffice or was undecided
    int found;                  // the states of the best controller
    int top;                    // search below it
    bool fewerFirst;
    bool unit;
    bool bisect = false;
    int step = 1;
    std::set<int> undecided;    // between refuted and top
};

#endif
