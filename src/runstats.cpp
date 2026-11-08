/**
 * Copyright Tom van Dijk
 */

#include <fstream>
#include <sstream>

#include <runstats.hpp>

namespace {

std::string quote(const std::string& s)
{
    std::string out = "\"";
    for (char c : s) {
        if (c == '"' || c == '\\') out += '\\';
        out += c;
    }
    return out + "\"";
}

const char* outcome(int result)
{
    switch (result) {
        case 1: return "suffice";
        case 0: return "do not suffice";
        case -2: return "too large";
        default: return "undecided";
    }
}

}  // namespace

bool RunStats::write(const std::string& path) const
{
    std::ostringstream o;
    o << "{\n  \"verdict\": " << (verdict.empty() ? "null" : quote(verdict)) << ",\n  \"times\": {";
    for (size_t i = 0; i < times.size(); i++) {
        o << (i ? ", " : "") << quote(times[i].first) << ": " << times[i].second;
    }
    const double total = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    o << (times.empty() ? "" : ", ") << "\"total\": " << total << "},\n  \"searches\": [";
    for (size_t i = 0; i < searches.size(); i++) {
        const MealyMinimizer::Search& s = searches[i];
        o << (i ? "," : "") << "\n    {\"moves\": " << quote(s.moves);
        if (!s.skipped.empty()) o << ", \"skipped\": " << quote(s.skipped);
        o << ", \"machine\": " << s.machine << ", \"quotient\": " << s.quotient << ", \"lower_bound\": "
          << s.lowerBound << ", \"found\": " << s.found << ", \"proven\": " << (s.proven ? "true" : "false")
          << ", \"calls\": [";
        for (size_t j = 0; j < s.calls.size(); j++) {
            const MealyMinimizer::Search::Call& call = s.calls[j];
            o << (j ? ", " : "") << "[" << call.states << ", " << quote(outcome(call.outcome)) << ", " << call.work << "]";
        }
        o << "]}";
    }
    o << (searches.empty() ? "" : "\n  ") << "],\n  \"candidates\": [";
    for (size_t i = 0; i < candidates.size(); i++) {
        const Candidate& c = candidates[i];
        o << (i ? "," : "") << "\n    {\"game\": " << quote(c.game) << ", \"encoding\": " << quote(c.encoding)
          << ", \"states_first\": " << (c.statesFirst ? "true" : "false") << ", \"all_candidates\": "
          << (c.allCandidates ? "true" : "false") << ", \"isop_candidate\": " << (c.isopCandidate ? "true" : "false")
          << ", \"reordered\": " << (c.reordered ? "true" : "false")
          << ", \"ands\": " << c.ands << ", \"latches\": " << c.latches << ", \"seconds\": " << c.seconds << "}";
    }
    o << (candidates.empty() ? "" : "\n  ") << "],\n  \"circuit\": ";
    if (ands < 0) {
        o << "null";
    } else {
        o << "{\"ands\": " << ands << ", \"latches\": " << latches << "}";
    }
    o << "\n}\n";
    std::ofstream file(path);
    file << o.str();
    return static_cast<bool>(file);
}
