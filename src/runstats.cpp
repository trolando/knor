/**
 * Copyright Tom van Dijk
 */

#include <cstdio>
#include <fstream>
#include <sstream>

#include <runstats.hpp>

namespace {

std::string quote(const std::string& s)
{
    std::string out = "\"";
    for (const char c : s) {
        if (c == '"' || c == '\\') {
            out += '\\';
            out += c;
        } else if (c == '\n') {
            out += "\\n";
        } else if (c == '\t') {
            out += "\\t";
        } else if (static_cast<unsigned char>(c) < 0x20) {
            char escaped[8];
            std::snprintf(escaped, sizeof(escaped), "\\u%04x", static_cast<unsigned char>(c));
            out += escaped;
        } else {
            out += c;
        }
    }
    return out + "\"";
}

const char* boolean(bool value) { return value ? "true" : "false"; }

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
    o << "{\n  \"schema\": " << SCHEMA << ",\n  \"configuration\": {";
    std::string stage;
    for (const Config::Setting& s : configuration) {
        if (s.stage != stage) {
            o << (stage.empty() ? "" : "}, ") << quote(s.stage) << ": {";
            stage = s.stage;
        } else {
            o << ", ";
        }
        o << quote(s.name) << ": ";
        switch (s.type) {
            case Config::Setting::Type::flag: o << boolean(s.value == "yes"); break;
            case Config::Setting::Type::number: o << s.value; break;
            case Config::Setting::Type::text: o << quote(s.value); break;
        }
    }
    o << (stage.empty() ? "" : "}") << "},\n  \"verdict\": " << (verdict.empty() ? "null" : quote(verdict))
      << ",\n  \"times\": {";
    for (size_t i = 0; i < times.size(); i++) {
        o << (i ? ", " : "") << quote(times[i].first) << ": " << times[i].second;
    }
    const double total = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    o << (times.empty() ? "" : ", ") << "\"total\": " << total << "},\n  \"extension\": ";
    if (!extension.ran) {
        o << "null";
    } else {
        const PermissiveStrategy::Stats& e = extension.stats;
        o << "{\"budget\": " << extension.budget << ", \"work\": " << e.work << ", \"stopped\": " << boolean(e.stopped)
          << ", \"vertices\": " << e.vertices << ", \"candidates\": " << e.candidates << ", \"added\": " << e.added << "}";
    }
    o << ",\n  \"searches\": [";
    for (size_t i = 0; i < searches.size(); i++) {
        const MealyMinimizer::Search& s = searches[i];
        o << (i ? "," : "") << "\n    {\"moves\": " << quote(s.moves);
        if (!s.skipped.empty()) o << ", \"skipped\": " << quote(s.skipped);
        if (!s.stop.empty()) o << ", \"stop\": " << quote(s.stop);
        o << ", \"machine\": " << s.machine << ", \"quotient\": " << s.quotient << ", \"lower_bound\": "
          << s.lowerBound << ", \"found\": " << s.found << ", \"proven\": " << boolean(s.proven) << ", \"calls\": [";
        for (size_t j = 0; j < s.calls.size(); j++) {
            const MealyMinimizer::Search::Call& call = s.calls[j];
            o << (j ? ", " : "") << "[" << call.states << ", " << quote(outcome(call.outcome)) << ", " << call.work << "]";
        }
        o << "]}";
    }
    o << (searches.empty() ? "" : "\n  ") << "],\n  \"candidates\": [";
    for (size_t i = 0; i < candidates.size(); i++) {
        const Candidate& c = candidates[i];
        o << (i ? "," : "") << "\n    {\"game\": " << quote(c.game) << ", \"kind\": " << quote(c.kind)
          << ", \"chosen\": " << boolean(c.chosen) << ", \"encoding\": " << quote(c.encoding)
          << ", \"states_first\": " << boolean(c.statesFirst) << ", \"all_candidates\": " << boolean(c.allCandidates)
          << ", \"isop_candidate\": " << boolean(c.isopCandidate) << ", \"reordered\": " << boolean(c.reorderWork > 0)
          << ", \"reorder_work\": " << c.reorderWork << ", \"reorder_used\": " << c.reorderUsed
          << ", \"ands\": " << c.ands << ", \"latches\": " << c.latches << ", \"seconds\": " << c.seconds << "}";
    }
    o << (candidates.empty() ? "" : "\n  ") << "],\n  \"selected\": ";
    if (selected < 0) {
        o << "null";
    } else {
        o << selected;
    }
    o << ",\n  \"circuit\": ";
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
