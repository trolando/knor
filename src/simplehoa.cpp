/**
 * Copyright Tom van Dijk
 */

/*
 * Parser for the extended HOA format (see simplehoa.h).
 *
 * The input is read into a single buffer that ends with a NUL sentinel, so
 * the hand-written lexer needs no bounds checks while scanning. A
 * recursive-descent parser builds the automaton in one pass. The body (the
 * bulk of a large automaton) is allocated in a bump arena: label trees are
 * hash-consed, so every atomic proposition, every negated proposition and
 * every recurring subexpression is stored once. Header data is staged in
 * standard containers and copied into malloc'd storage when the automaton
 * is complete, so a failed or aborted parse never touches the output.
 */

#include <algorithm>
#include <array>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <new>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "hoainternal.hpp"
#include "simplehoa.h"

using namespace hoa;

namespace {

// ---------------------------------------------------------------------------
// Hash-consed label trees
// ---------------------------------------------------------------------------

class LabelFactory final {
public:
    explicit LabelFactory(Arena& storage) : arena(storage), table(INITIAL_TABLE, nullptr)
    {
        constants[0] = node(NT_BOOL, nullptr, nullptr, 0, nullptr);
        constants[1] = node(NT_BOOL, nullptr, nullptr, 1, nullptr);
    }

    BTree* constant(bool value) { return constants[value ? 1 : 0]; }

    BTree* ap(int id)
    {
        if (static_cast<size_t>(id) >= aps.size()) aps.resize(id + 1, nullptr);
        BTree*& n = aps[id];
        if (n == nullptr) n = node(NT_AP, nullptr, nullptr, id, nullptr);
        return n;
    }

    BTree* alias(char* name) { return node(NT_ALIAS, nullptr, nullptr, -1, name); }

    BTree* negate(BTree* child)
    {
        if (child->type != NT_AP) return compound(NT_NOT, child, nullptr);
        if (static_cast<size_t>(child->id) >= negatedAps.size()) negatedAps.resize(child->id + 1, nullptr);
        BTree*& n = negatedAps[child->id];
        if (n == nullptr) n = node(NT_NOT, child, nullptr, -1, nullptr);
        return n;
    }

    BTree* binary(NodeType type, BTree* left, BTree* right) { return compound(type, left, right); }

private:
    static constexpr size_t INITIAL_TABLE = size_t(1) << 12;

    BTree* node(NodeType type, BTree* left, BTree* right, int id, char* alias)
    {
        BTree* n = arena.allocate<BTree>(1);
        n->left = left;
        n->right = right;
        n->alias = alias;
        n->id = id;
        n->type = type;
        return n;
    }

    static size_t hash(NodeType type, const BTree* left, const BTree* right)
    {
        uint64_t h = reinterpret_cast<uintptr_t>(left);
        h ^= reinterpret_cast<uintptr_t>(right) * 0xC2B2AE3D27D4EB4FULL;
        h ^= static_cast<uint64_t>(type) << 1;
        h *= 0x9E3779B97F4A7C15ULL;
        return static_cast<size_t>(h ^ (h >> 32));
    }

    // Open-addressing hash table of all operator nodes (linear probing).
    BTree* compound(NodeType type, BTree* left, BTree* right)
    {
        size_t mask = table.size() - 1;
        size_t i = hash(type, left, right) & mask;
        while (BTree* n = table[i]) {
            if (n->left == left && n->right == right && n->type == type) return n;
            i = (i + 1) & mask;
        }
        BTree* n = node(type, left, right, -1, nullptr);
        table[i] = n;
        if (++count * 2 > table.size()) rehash();
        return n;
    }

    void rehash()
    {
        std::vector<BTree*> old(table.size() * 2, nullptr);
        old.swap(table);
        const size_t mask = table.size() - 1;
        for (BTree* n : old) {
            if (n == nullptr) continue;
            size_t i = hash(n->type, n->left, n->right) & mask;
            while (table[i] != nullptr) i = (i + 1) & mask;
            table[i] = n;
        }
    }

    Arena& arena;
    BTree* constants[2];
    std::vector<BTree*> aps;
    std::vector<BTree*> negatedAps;
    std::vector<BTree*> table;
    size_t count = 0;
};

// Labels in the body of large automata are highly repetitive, so the parser
// maps the text of every label it has parsed to the resulting tree, and looks
// up the text of a label before parsing it.
class LabelTexts final {
public:
    LabelTexts() : table(INITIAL_TABLE) {}

    static uint64_t hash(const char* s, size_t n)
    {
        uint64_t h = 0x9E3779B97F4A7C15ULL ^ n;
        if (n >= 8) {
            const char* last = s + n - 8;
            for (; s < last; s += 8) h = mix(h, load8(s));
            return mix(h, load8(last));     // the final 8 bytes, possibly overlapping
        }
        uint64_t w = 0;
        for (size_t i = 0; i < n; i++) w |= static_cast<uint64_t>(static_cast<unsigned char>(s[i])) << (8 * i);
        return mix(h, w);
    }

    BTree* find(const char* text, size_t length, uint64_t h) const
    {
        const size_t mask = table.size() - 1;
        for (size_t i = h & mask; table[i].root != nullptr; i = (i + 1) & mask) {
            const Entry& e = table[i];
            if (e.hash == h && e.length == length && memcmp(e.text, text, length) == 0) return e.root;
        }
        return nullptr;
    }

    void insert(const char* text, size_t length, uint64_t h, BTree* root)
    {
        place(Entry{text, length, h, root});
        if (++count * 2 > table.size()) {
            std::vector<Entry> old(table.size() * 2);
            old.swap(table);
            for (const Entry& e : old) {
                if (e.root != nullptr) place(e);
            }
        }
    }

private:
    static constexpr size_t INITIAL_TABLE = size_t(1) << 10;

    static uint64_t load8(const char* s)
    {
        uint64_t w;
        memcpy(&w, s, 8);
        return w;
    }

    static uint64_t mix(uint64_t h, uint64_t w)
    {
        h = (h ^ w) * 0xFF51AFD7ED558CCDULL;
        return h ^ (h >> 32);
    }

    struct Entry {
        const char* text = nullptr;
        size_t length = 0;
        uint64_t hash = 0;
        BTree* root = nullptr;
    };

    void place(const Entry& e)
    {
        const size_t mask = table.size() - 1;
        size_t i = e.hash & mask;
        while (table[i].root != nullptr) i = (i + 1) & mask;
        table[i] = e;
    }

    std::vector<Entry> table;
    size_t count = 0;
};

// ---------------------------------------------------------------------------
// Lexer
// ---------------------------------------------------------------------------

struct ParseError {
    const char* position;
    std::string message;
};

// Thrown when the current automaton is ended by --ABORT--.
struct Aborted {};

enum class Tok : uint8_t {
    End,        // end of input
    Int,        // value
    Bool,       // value 1 (t) or 0 (f)
    Ident,      // text/length
    String,     // text/length (raw, between the quotes), escaped if it contains '\'
    AName,      // text/length (without the '@')
    Header,     // text/length (without the ':')
    LParen, RParen, LBrace, RBrace, LBracket, RBracket,
    And, Or, Not,
    Body,       // --BODY--
    EndMark,    // --END--
};

enum : uint8_t {
    CC_SPACE = 1,
    CC_IDSTART = 2,     // first character of an identifier
    CC_IDCHAR = 4,      // other characters of identifiers and alias names
};

constexpr std::array<uint8_t, 256> makeCharClasses()
{
    std::array<uint8_t, 256> cc{};
    cc[' '] = cc['\t'] = cc['\n'] = cc['\r'] = CC_SPACE;
    for (int c = 0; c < 256; c++) {
        const bool alpha = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
        const bool digit = c >= '0' && c <= '9';
        if (alpha || c == '_') cc[c] |= CC_IDSTART;
        // HOA v1 does not allow '.'; like Spot, accept it, as HOA v1.1 version
        // strings and header names (spot.highlight.edges:) contain it.
        if (alpha || digit || c == '_' || c == '-' || c == '.') cc[c] |= CC_IDCHAR;
    }
    return cc;
}

constexpr std::array<uint8_t, 256> charClasses = makeCharClasses();

constexpr uint8_t NOT_SINGLE = 0xff;

constexpr std::array<uint8_t, 256> makeSingleCharTokens()
{
    std::array<uint8_t, 256> t{};
    for (auto& k : t) k = NOT_SINGLE;
    t['['] = static_cast<uint8_t>(Tok::LBracket);
    t[']'] = static_cast<uint8_t>(Tok::RBracket);
    t['&'] = static_cast<uint8_t>(Tok::And);
    t['|'] = static_cast<uint8_t>(Tok::Or);
    t['!'] = static_cast<uint8_t>(Tok::Not);
    t['('] = static_cast<uint8_t>(Tok::LParen);
    t[')'] = static_cast<uint8_t>(Tok::RParen);
    t['{'] = static_cast<uint8_t>(Tok::LBrace);
    t['}'] = static_cast<uint8_t>(Tok::RBrace);
    return t;
}

constexpr std::array<uint8_t, 256> singleCharTokens = makeSingleCharTokens();

inline bool hasClass(char c, uint8_t cls)
{
    return (charClasses[static_cast<unsigned char>(c)] & cls) != 0;
}

class Lexer final {
public:
    // Lex first..last; the character at last must be a NUL sentinel.
    Lexer(const char* first, const char* last) : begin(first), end(last), p(first) {}

    Tok kind = Tok::End;
    int value = 0;
    const char* text = nullptr;
    size_t length = 0;
    bool escaped = false;
    const char* start = nullptr;    // first character of the current token

    // Advance to the next token. The common tokens of edges (punctuation and
    // small integers) are handled inline; everything else out of line.
    void next()
    {
        while (hasClass(*p, CC_SPACE)) p++;
        start = p;
        const unsigned char c = static_cast<unsigned char>(*p);
        const uint8_t single = singleCharTokens[c];
        if (single != NOT_SINGLE) {
            p++;
            kind = static_cast<Tok>(single);
            return;
        }
        if (c >= '0' && c <= '9') {
            p++;
            kind = Tok::Int;
            value = c - '0';
            if (c != '0' && *p >= '0' && *p <= '9') lexIntTail();
            return;
        }
        nextSlow();
    }

    bool isHeader(const char* name) const
    {
        return kind == Tok::Header && length == strlen(name) && memcmp(text, name, length) == 0;
    }

    bool isIdent(const char* name) const
    {
        return kind == Tok::Ident && length == strlen(name) && memcmp(text, name, length) == 0;
    }

    std::string tokenText() const { return std::string(text, length); }

    // The input after the current token.
    const char* rest() const { return p; }
    size_t restLength() const { return static_cast<size_t>(end - p); }

    // Continue lexing at the given position of the input.
    void resume(const char* at)
    {
        p = at;
        next();
    }

    // The current string token with its escape sequences resolved.
    std::string stringValue() const
    {
        if (!escaped) return tokenText();
        std::string s;
        s.reserve(length);
        for (size_t i = 0; i < length; i++) {
            if (text[i] == '\\') i++;
            s.push_back(text[i]);
        }
        return s;
    }

    [[noreturn]] static void error(const char* at, std::string message)
    {
        throw ParseError{at, std::move(message)};
    }

    // 1-based line and column of a position in the input.
    std::pair<long, long> lineColumn(const char* at) const
    {
        long line = 1;
        const char* lineStart = begin;
        for (const char* q = begin; q < at; q++) {
            if (*q == '\n') {
                line++;
                lineStart = q + 1;
            }
        }
        return {line, static_cast<long>(at - lineStart) + 1};
    }

private:
    [[gnu::noinline]] void nextSlow()
    {
        switch (*p) {
            case '/':
                // skip all comments here, so the next() below lands on a token
                do {
                    if (p[1] != '*') error(p, "unexpected character '/'");
                    skipComment();
                    while (hasClass(*p, CC_SPACE)) p++;
                } while (*p == '/');
                next();
                return;
            case '"': lexString(); return;
            case '@': lexAName(); return;
            case '-': lexMarker(); return;
            case '\0':
                if (p != end) error(p, "unexpected NUL character");
                kind = Tok::End;
                return;
            default:
                if (!hasClass(*p, CC_IDSTART)) error(p, std::string("unexpected character '") + *p + "'");
                lexIdent();
                return;
        }
    }

    // Comments nest: /* a /* b */ c */ is one comment.
    void skipComment()
    {
        const char* open = p;
        int depth = 1;
        p += 2;
        while (depth > 0) {
            if (p[0] == '*' && p[1] == '/') {
                depth--;
                p += 2;
            } else if (p[0] == '/' && p[1] == '*') {
                depth++;
                p += 2;
            } else if (p[0] == '\0' && p == end) {
                error(open, "unterminated comment");
            } else {
                p++;
            }
        }
    }

    // The remaining digits of an integer whose first digit is in value.
    [[gnu::noinline]] void lexIntTail()
    {
        int64_t v = value;
        while (*p >= '0' && *p <= '9') {
            v = v * 10 + (*p++ - '0');
            if (v > INT_MAX) error(start, "integer too large");
        }
        value = static_cast<int>(v);
    }

    void lexIdent()
    {
        const char* s = p++;
        while (hasClass(*p, CC_IDCHAR)) p++;
        text = s;
        length = static_cast<size_t>(p - s);
        if (*p == ':') {
            p++;
            kind = Tok::Header;
        } else if (length == 1 && (*s == 't' || *s == 'f')) {
            kind = Tok::Bool;
            value = *s == 't' ? 1 : 0;
        } else {
            kind = Tok::Ident;
        }
    }

    void lexAName()
    {
        const char* s = ++p;
        while (hasClass(*p, CC_IDCHAR)) p++;
        if (p == s) error(start, "expected an alias name after '@'");
        kind = Tok::AName;
        text = s;
        length = static_cast<size_t>(p - s);
    }

    void lexString()
    {
        const char* s = ++p;
        escaped = false;
        for (;;) {
            const char c = *p;
            if (c == '"') break;
            if (c == '\\') {
                escaped = true;
                p++;
            }
            if (*p == '\0' && p == end) error(start, "unterminated string");
            p++;
        }
        kind = Tok::String;
        text = s;
        length = static_cast<size_t>(p - s);
        p++;
    }

    bool matchMarker(const char* marker, size_t n)
    {
        if (strncmp(p, marker, n) != 0) return false;
        p += n;
        return true;
    }

    void lexMarker()
    {
        if (matchMarker("--BODY--", 8)) {
            kind = Tok::Body;
        } else if (matchMarker("--END--", 7)) {
            kind = Tok::EndMark;
        } else if (matchMarker("--ABORT--", 9)) {
            throw Aborted{};
        } else {
            error(p, "unexpected character '-'");
        }
    }

    const char* begin;
    const char* end;
    const char* p;
};

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

// Nesting limit for parentheses and negations, so that hostile input cannot
// exhaust the (Lace worker) stack.
constexpr int MAX_NESTING = 1000;

class AutomatonParser final {
public:
    explicit AutomatonParser(Lexer& lexer) : lex(lexer), arena(new Arena), labels(*arena) {}

    AutomatonParser(const AutomatonParser&) = delete;
    AutomatonParser& operator=(const AutomatonParser&) = delete;

    // The acceptance signatures of edges are malloc'd; free them unless they
    // were handed over by store().
    ~AutomatonParser()
    {
        if (stored) return;
        for (const StagedState& s : states) {
            for (int j = 0; j < s.noTrans; j++) free(s.transitions[j].accSig);
        }
        for (const Transition& t : edges) free(t.accSig);
    }

    // Parse one automaton; the first token must be current.
    void parse()
    {
        if (!lex.isHeader("HOA")) lex.error(lex.start, "expected \"HOA:\" at the start of the automaton");
        lex.next();
        if (lex.kind != Tok::Ident) lex.error(lex.start, "expected a format version after \"HOA:\"");
        version = lex.tokenText();
        lex.next();
        while (lex.kind == Tok::Header) parseHeaderItem();
        if (lex.kind != Tok::Body) lex.error(lex.start, "expected a header item or \"--BODY--\"");
        startBody();
        lex.next();
        while (lex.isHeader("State")) parseState();
        if (lex.kind != Tok::EndMark) lex.error(lex.start, "expected \"State:\" or \"--END--\"");
        finishBody();
        // Only the end of the input or another automaton (which is ignored)
        // may follow.
        try {
            lex.next();
        } catch (const Aborted&) {
            return;
        }
        if (lex.kind != Tok::End && !lex.isHeader("HOA")) lex.error(lex.start, "unexpected input after \"--END--\"");
    }

    // Move the parsed automaton into data (which is overwritten).
    void store(HoaData* data);

private:
    enum HeaderItem { H_STATES, H_AP, H_CNTAP, H_ACCEPTANCE, H_ACCNAME, H_TOOL, H_NAME, H_COUNT };

    void once(HeaderItem item, const char* header)
    {
        if (seen[item]) lex.error(lex.start, std::string("more than one \"") + header + ":\" header item");
        seen[item] = true;
    }

    int expectInt(const char* what)
    {
        if (lex.kind != Tok::Int) lex.error(lex.start, std::string("expected ") + what);
        const int v = lex.value;
        lex.next();
        return v;
    }

    void expect(Tok kind, const char* what)
    {
        if (lex.kind != kind) lex.error(lex.start, std::string("expected ") + what);
        lex.next();
    }

    void parseHeaderItem()
    {
        const char* item = lex.start;
        if (lex.isHeader("States")) {
            once(H_STATES, "States");
            lex.next();
            declaredStates = expectInt("the number of states");
        } else if (lex.isHeader("Start")) {
            lex.next();
            parseStateConj(start);
        } else if (lex.isHeader("AP")) {
            once(H_AP, "AP");
            lex.next();
            const int n = expectInt("the number of atomic propositions");
            while (lex.kind == Tok::String) {
                aps.push_back(lex.stringValue());
                lex.next();
            }
            if (aps.size() != static_cast<size_t>(n)) {
                lex.error(item, "\"AP:\" declares " + std::to_string(n) + " atomic propositions but lists "
                                + std::to_string(aps.size()));
            }
            apBound = n;
        } else if (lex.isHeader("controllable-AP")) {
            once(H_CNTAP, "controllable-AP");
            cntAPsPos = item;
            lex.next();
            while (lex.kind == Tok::Int) {
                cntAPs.push_back(lex.value);
                lex.next();
            }
        } else if (lex.isHeader("Alias")) {
            lex.next();
            parseAlias();
        } else if (lex.isHeader("Acceptance")) {
            once(H_ACCEPTANCE, "Acceptance");
            lex.next();
            noAccSets = expectInt("the number of acceptance sets");
            acc = parseAccOr(0);
        } else if (lex.isHeader("acc-name")) {
            once(H_ACCNAME, "acc-name");
            lex.next();
            if (lex.kind != Tok::Ident) lex.error(lex.start, "expected an acceptance name");
            accName = lex.tokenText();
            hasAccName = true;
            lex.next();
            for (;;) {
                if (lex.kind == Tok::Ident) accNameParameters.push_back(lex.tokenText());
                else if (lex.kind == Tok::Int) accNameParameters.push_back(std::to_string(lex.value));
                else if (lex.kind == Tok::Bool) accNameParameters.push_back(lex.value ? "True" : "False");
                else break;
                lex.next();
            }
        } else if (lex.isHeader("tool")) {
            once(H_TOOL, "tool");
            lex.next();
            if (lex.kind != Tok::String) lex.error(lex.start, "expected the tool name");
            tool = lex.stringValue();
            hasTool = true;
            lex.next();
            if (lex.kind == Tok::String) {
                toolVersion = lex.stringValue();
                hasToolVersion = true;
                lex.next();
            }
        } else if (lex.isHeader("name")) {
            once(H_NAME, "name");
            lex.next();
            if (lex.kind != Tok::String) lex.error(lex.start, "expected the automaton name");
            name = lex.stringValue();
            hasName = true;
            lex.next();
        } else if (lex.isHeader("properties")) {
            lex.next();
            for (;;) {
                if (lex.kind == Tok::Ident) {
                    properties.push_back(lex.tokenText());
                } else if (lex.kind == Tok::Not) {
                    // HOA v1.1: a property that is known not to hold
                    lex.next();
                    if (lex.kind != Tok::Ident) lex.error(lex.start, "expected a property after '!'");
                    properties.push_back("!" + lex.tokenText());
                } else {
                    break;
                }
                lex.next();
            }
        } else if (lex.isHeader("HOA") || lex.isHeader("State")) {
            lex.error(item, "unexpected \"" + lex.tokenText() + ":\" in the header");
        } else {
            // Unknown header item: skip its arguments.
            std::string header = lex.tokenText();
            if (header[0] >= 'A' && header[0] <= 'Z') {
                auto lc = lex.lineColumn(item);
                fprintf(stderr, "Warning: ignoring unsupported header item \"%s:\" [line %ld, column %ld]\n",
                        header.c_str(), lc.first, lc.second);
            }
            lex.next();
            while (lex.kind == Tok::Int || lex.kind == Tok::Bool || lex.kind == Tok::String
                   || lex.kind == Tok::Ident) {
                lex.next();
            }
        }
    }

    void parseAlias()
    {
        if (lex.kind != Tok::AName) lex.error(lex.start, "expected an alias name");
        const char* at = lex.start;
        std::string aliasName = lex.tokenText();
        if (aliasIndex.count(aliasName) != 0) lex.error(at, "alias @" + aliasName + " is defined twice");
        lex.next();
        BTree* expr = parseLabelOr(0);
        char* shared = arena->copyString(aliasName);
        aliasIndex.emplace(aliasName, aliases.size());
        aliases.push_back(AliasEntry{std::move(aliasName), expr, labels.alias(shared)});
    }

    // state-conj: INT ("&" INT)*
    void parseStateConj(std::vector<int>& out)
    {
        out.push_back(expectState());
        while (lex.kind == Tok::And) {
            lex.next();
            out.push_back(expectState());
        }
    }

    int expectState()
    {
        if (lex.kind != Tok::Int) lex.error(lex.start, "expected a state number");
        const int s = lex.value;
        if (inBody && declaredStates >= 0 && s >= declaredStates) {
            lex.error(lex.start, "state " + std::to_string(s) + " is out of range (States: "
                                 + std::to_string(declaredStates) + ")");
        }
        if (s > maxState) {
            maxState = s;
            maxStatePos = lex.start;
        }
        lex.next();
        return s;
    }

    // acc-sig: "{" INT* "}"; the current token is "{".
    void parseAccSig(std::vector<int>& out)
    {
        lex.next();
        while (lex.kind == Tok::Int) {
            if (lex.value >= noAccSets) {
                lex.error(lex.start, "acceptance set " + std::to_string(lex.value) + " is out of range (Acceptance: "
                                     + std::to_string(noAccSets) + ")");
            }
            out.push_back(lex.value);
            lex.next();
        }
        expect(Tok::RBrace, "an acceptance set or '}'");
    }

    // label-expr, with & binding stronger than |, both left-associative.
    BTree* parseLabelOr(int depth)
    {
        BTree* left = parseLabelAnd(depth);
        while (lex.kind == Tok::Or) {
            lex.next();
            BTree* right = parseLabelAnd(depth);
            left = labels.binary(NT_OR, left, right);
        }
        return left;
    }

    BTree* parseLabelAnd(int depth)
    {
        BTree* left = parseLabelAtom(depth);
        while (lex.kind == Tok::And) {
            lex.next();
            BTree* right = parseLabelAtom(depth);
            left = labels.binary(NT_AND, left, right);
        }
        return left;
    }

    BTree* parseLabelAtom(int depth)
    {
        switch (lex.kind) {
            case Tok::Int: {
                const int id = lex.value;
                if (id >= apBound) {
                    lex.error(lex.start, "atomic proposition " + std::to_string(id) + " is out of range (AP: "
                                         + std::to_string(aps.size()) + ")");
                }
                if (!inBody && id > maxHeaderAP) {
                    maxHeaderAP = id;
                    maxHeaderAPPos = lex.start;
                }
                lex.next();
                return labels.ap(id);
            }
            case Tok::Not: {
                if (depth >= MAX_NESTING) lex.error(lex.start, "label expression nested too deeply");
                lex.next();
                return labels.negate(parseLabelAtom(depth + 1));
            }
            case Tok::Bool: {
                const bool value = lex.value != 0;
                lex.next();
                return labels.constant(value);
            }
            case Tok::AName: {
                auto it = aliasIndex.find(lex.tokenText());
                if (it == aliasIndex.end()) lex.error(lex.start, "alias @" + lex.tokenText() + " is not defined");
                lex.next();
                return aliases[it->second].node;
            }
            case Tok::LParen: {
                if (depth >= MAX_NESTING) lex.error(lex.start, "label expression nested too deeply");
                lex.next();
                BTree* e = parseLabelOr(depth + 1);
                expect(Tok::RParen, "')'");
                return e;
            }
            default:
                lex.error(lex.start, "expected a label expression");
        }
    }

    // acceptance-cond, with & binding stronger than |, both left-associative.
    TreePtr parseAccOr(int depth)
    {
        TreePtr left = parseAccAnd(depth);
        while (lex.kind == Tok::Or) {
            lex.next();
            TreePtr right = parseAccAnd(depth);
            left = mallocNode(NT_OR, std::move(left), std::move(right), -1);
        }
        return left;
    }

    TreePtr parseAccAnd(int depth)
    {
        TreePtr left = parseAccAtom(depth);
        while (lex.kind == Tok::And) {
            lex.next();
            TreePtr right = parseAccAtom(depth);
            left = mallocNode(NT_AND, std::move(left), std::move(right), -1);
        }
        return left;
    }

    TreePtr parseAccAtom(int depth)
    {
        if (lex.kind == Tok::Bool) {
            const int value = lex.value;
            lex.next();
            return mallocNode(NT_BOOL, nullptr, nullptr, value);
        }
        if (lex.kind == Tok::LParen) {
            if (depth >= MAX_NESTING) lex.error(lex.start, "acceptance condition nested too deeply");
            lex.next();
            TreePtr e = parseAccOr(depth + 1);
            expect(Tok::RParen, "')'");
            return e;
        }
        NodeType type;
        if (lex.isIdent("Fin")) type = NT_FIN;
        else if (lex.isIdent("Inf")) type = NT_INF;
        else lex.error(lex.start, "expected Fin, Inf, t, f or '(' in the acceptance condition");
        lex.next();
        expect(Tok::LParen, "'('");
        const bool negated = lex.kind == Tok::Not;
        if (negated) lex.next();
        if (lex.kind == Tok::Int && lex.value >= noAccSets) {
            lex.error(lex.start, "acceptance set " + std::to_string(lex.value) + " is out of range (Acceptance: "
                                 + std::to_string(noAccSets) + ")");
        }
        TreePtr set = mallocNode(NT_SET, nullptr, nullptr, expectInt("an acceptance set"));
        if (negated) set = mallocNode(NT_NOT, std::move(set), nullptr, -1);
        expect(Tok::RParen, "')'");
        return mallocNode(type, std::move(set), nullptr, -1);
    }

    // Header checks that need the complete header.
    void startBody()
    {
        if (!seen[H_ACCEPTANCE]) lex.error(lex.start, "missing \"Acceptance:\" header item");
        apBound = static_cast<int>(aps.size());
        if (maxHeaderAP >= apBound) {
            lex.error(maxHeaderAPPos, "atomic proposition " + std::to_string(maxHeaderAP) + " is out of range (AP: "
                                     + std::to_string(aps.size()) + ")");
        }
        for (int c : cntAPs) {
            if (c >= apBound) {
                lex.error(cntAPsPos, "controllable atomic proposition " + std::to_string(c)
                                     + " is out of range (AP: " + std::to_string(aps.size()) + ")");
            }
        }
        if (declaredStates >= 0 && maxState >= declaredStates) {
            lex.error(maxStatePos, "state " + std::to_string(maxState) + " is out of range (States: "
                                   + std::to_string(declaredStates) + ")");
        }
        inBody = true;
        if (declaredStates >= 0) states.resize(declaredStates);
    }

    void parseState()
    {
        lex.next();
        BTree* stateLabel = nullptr;
        if (lex.kind == Tok::LBracket) stateLabel = parseBodyLabel();
        const char* idPos = lex.start;
        const int id = expectState();
        if (static_cast<size_t>(id) >= states.size()) states.resize(static_cast<size_t>(id) + 1);
        StagedState& s = states[id];
        if (s.defined) lex.error(idPos, "state " + std::to_string(id) + " is defined twice");
        s.defined = true;
        s.label = stateLabel;
        if (lex.kind == Tok::String) {
            s.name = arena->copyString(lex.stringValue());
            lex.next();
        }
        if (lex.kind == Tok::LBrace) {
            ints.clear();
            parseAccSig(ints);
            s.accSig = ints;
        }

        // edges: label? state-conj acc-sig?
        edges.clear();
        size_t unlabelled = 0;
        while (lex.kind == Tok::LBracket || lex.kind == Tok::Int) {
            Transition t;
            t.accSig = nullptr;
            t.noAccSig = 0;
            if (lex.kind == Tok::LBracket) {
                if (stateLabel != nullptr) lex.error(lex.start, "edge label in a state that has a label");
                t.label = parseBodyLabel();
            } else {
                t.label = nullptr;
                unlabelled++;
            }
            if (lex.kind == Tok::Int) {
                // fast path: a single successor
                const int succ = expectState();
                if (lex.kind == Tok::And) {
                    ints.clear();
                    ints.push_back(succ);
                    while (lex.kind == Tok::And) {
                        lex.next();
                        ints.push_back(expectState());
                    }
                    t.successors = arena->copy(ints.data(), ints.size());
                    t.noSucc = static_cast<int>(ints.size());
                } else {
                    t.successors = arena->allocate<int>(1);
                    t.successors[0] = succ;
                    t.noSucc = 1;
                }
            } else {
                lex.error(lex.start, "expected a state number");
            }
            if (lex.kind == Tok::LBrace) {
                ints.clear();
                parseAccSig(ints);
                edges.push_back(t);
                edges.back().accSig = mallocArray(ints);
                edges.back().noAccSig = static_cast<int>(ints.size());
            } else {
                edges.push_back(t);
            }
        }

        if (unlabelled != 0 && stateLabel == nullptr) {
            if (unlabelled != edges.size()) {
                lex.error(idPos, "state " + std::to_string(id) + " mixes labelled and unlabelled edges");
            }
            labelImplicitEdges(id, idPos);
        }
        s.transitions = arena->copy(edges.data(), edges.size());
        s.noTrans = static_cast<int>(edges.size());
        edges.clear();
    }

    // A label in the body; the current token is '['. Its text up to the
    // first ']' is looked up first: if it equals the text of a label that
    // was parsed before, it is that label (a ']' inside a comment yields a
    // text with an unclosed comment, which is never in the table).
    // Otherwise, a label in disjunctive normal form is parsed directly from
    // the text, and any other label by the general parser.
    BTree* parseBodyLabel()
    {
        const char* text = lex.rest();
        const char* close = static_cast<const char*>(memchr(text, ']', lex.restLength()));
        if (close == nullptr) {
            lex.next();
            BTree* label = parseLabelOr(0);
            expect(Tok::RBracket, "']'");
            return label;
        }
        const size_t length = static_cast<size_t>(close - text);
        const uint64_t h = LabelTexts::hash(text, length);
        if (BTree* known = labelTexts.find(text, length, h)) {
            lex.resume(close + 1);
            return known;
        }
        BTree* label = parseDnf(text, close);
        if (label != nullptr) {
            lex.resume(close + 1);
        } else {
            lex.next();
            label = parseLabelOr(0);
            if (lex.kind != Tok::RBracket || lex.start != close) {
                expect(Tok::RBracket, "']'");
                return label;
            }
            lex.next();
        }
        labelTexts.insert(text, length, h, label);
        return label;
    }

    // Parse text..end as a disjunction of cubes, each cube looked up by its
    // text first. Returns nullptr, without side effects other than creating
    // label nodes, if the text has any other form (or is invalid); the
    // general parser then handles it, so both give the same tree.
    BTree* parseDnf(const char* text, const char* end)
    {
        BTree* result = nullptr;
        for (const char* q = text;;) {
            const char* bar = static_cast<const char*>(memchr(q, '|', static_cast<size_t>(end - q)));
            const char* cubeEnd = bar != nullptr ? bar : end;
            const size_t length = static_cast<size_t>(cubeEnd - q);
            const uint64_t h = LabelTexts::hash(q, length);
            BTree* cube = cubeTexts.find(q, length, h);
            if (cube == nullptr) {
                cube = parseCube(q, cubeEnd);
                if (cube == nullptr) return nullptr;
                cubeTexts.insert(q, length, h, cube);
            }
            result = result == nullptr ? cube : labels.binary(NT_OR, result, cube);
            if (bar == nullptr) return result;
            q = bar + 1;
        }
    }

    // Parse text..end as a conjunction of literals: [!] (INT | t | f).
    BTree* parseCube(const char* q, const char* end)
    {
        BTree* result = nullptr;
        for (;;) {
            while (q != end && hasClass(*q, CC_SPACE)) q++;
            bool negated = false;
            if (q != end && *q == '!') {
                negated = true;
                q++;
                while (q != end && hasClass(*q, CC_SPACE)) q++;
            }
            if (q == end) return nullptr;
            BTree* literal;
            if (*q >= '1' && *q <= '9') {
                int64_t id = *q++ - '0';
                while (q != end && *q >= '0' && *q <= '9') {
                    id = id * 10 + (*q++ - '0');
                    if (id >= apBound) return nullptr;
                }
                if (id >= apBound) return nullptr;
                literal = labels.ap(static_cast<int>(id));
            } else if (*q == '0') {
                if (apBound == 0) return nullptr;
                q++;
                literal = labels.ap(0);
            } else if (*q == 't' || *q == 'f') {
                literal = labels.constant(*q == 't');
                q++;
            } else {
                return nullptr;
            }
            // the literal must end here: "01", "0a", "tt", "t:" are not cubes
            if (q != end && (hasClass(*q, CC_IDCHAR) || *q == ':')) return nullptr;
            if (negated) literal = labels.negate(literal);
            result = result == nullptr ? literal : labels.binary(NT_AND, result, literal);
            while (q != end && hasClass(*q, CC_SPACE)) q++;
            if (q == end) return result;
            if (*q != '&') return nullptr;
            q++;
        }
    }

    // With implicit labels, edge i of a state carries the cube of valuation i.
    void labelImplicitEdges(int id, const char* idPos)
    {
        const int n = apBound;
        if (n >= 31 || edges.size() != (size_t(1) << n)) {
            lex.error(idPos, "state " + std::to_string(id) + " has " + std::to_string(edges.size())
                             + " unlabelled edges, but implicit labels require 2^" + std::to_string(n));
        }
        for (size_t i = 0; i < edges.size(); i++) {
            BTree* cube = nullptr;
            for (int j = 0; j < n; j++) {
                BTree* lit = labels.ap(j);
                if ((i & (size_t(1) << j)) == 0) lit = labels.negate(lit);
                cube = cube == nullptr ? lit : labels.binary(NT_AND, cube, lit);
            }
            edges[i].label = cube == nullptr ? labels.constant(true) : cube;
        }
    }

    // Every state that is used (initial or successor) must have a "State:"
    // entry; other states without one are kept as dead ends (as Spot does).
    void finishBody()
    {
        const size_t n = declaredStates >= 0 ? static_cast<size_t>(declaredStates) : static_cast<size_t>(maxState + 1);
        states.resize(n);
        auto mustBeDefined = [&](int s) {
            if (!states[s].defined) {
                lex.error(lex.start, "state " + std::to_string(s) + " is used but has no \"State:\" entry");
            }
        };
        for (int s : start) mustBeDefined(s);
        for (const StagedState& s : states) {
            for (int j = 0; j < s.noTrans; j++) {
                const Transition& t = s.transitions[j];
                for (int k = 0; k < t.noSucc; k++) mustBeDefined(t.successors[k]);
            }
        }
        size_t undefined = 0;
        for (const StagedState& s : states) {
            if (!s.defined) undefined++;
        }
        if (undefined != 0) {
            auto lc = lex.lineColumn(lex.start);
            fprintf(stderr, "Warning: %zu unused state(s) without a \"State:\" entry [line %ld, column %ld]\n",
                    undefined, lc.first, lc.second);
        }
    }

    struct StagedState {
        bool defined = false;
        char* name = nullptr;
        BTree* label = nullptr;
        std::vector<int> accSig;
        Transition* transitions = nullptr;
        int noTrans = 0;
    };

    struct AliasEntry {
        std::string name;
        BTree* expr;
        BTree* node;        // the NT_ALIAS node shared by all uses
    };

    Lexer& lex;
    std::unique_ptr<Arena> arena;
    LabelFactory labels;

    bool seen[H_COUNT] = {};
    std::string version, accName, tool, toolVersion, name;
    bool hasAccName = false, hasTool = false, hasToolVersion = false, hasName = false;
    std::vector<std::string> aps, accNameParameters, properties;
    std::vector<int> start, cntAPs;
    const char* cntAPsPos = nullptr;
    std::vector<AliasEntry> aliases;
    std::unordered_map<std::string, size_t> aliasIndex;
    int noAccSets = -1;
    TreePtr acc;
    int declaredStates = -1;

    // Atomic propositions are range-checked as soon as "AP:" is known; those
    // in aliases that precede "AP:" are checked at --BODY--.
    int apBound = INT_MAX;
    int maxHeaderAP = -1;
    const char* maxHeaderAPPos = nullptr;
    int maxState = -1;
    const char* maxStatePos = nullptr;
    bool inBody = false;

    std::vector<StagedState> states;
    std::vector<Transition> edges;  // edges of the current state
    std::vector<int> ints;          // scratch list
    LabelTexts labelTexts;
    LabelTexts cubeTexts;
    bool stored = false;
};

void AutomatonParser::store(HoaData* data)
{
    // First allocate all malloc'd parts. If that fails, out holds only
    // malloc'd memory, which resetHoa frees, and the arena and the edges
    // remain with the parser.
    HoaData out;
    defaultsHoa(&out);
    try {
        out.version = mallocString(version);
        if (hasName) out.name = mallocString(name);
        if (hasAccName) out.accNameID = mallocString(accName);
        if (hasTool) out.toolName = mallocString(tool);
        if (hasToolVersion) out.toolVersion = mallocString(toolVersion);

        auto strings = [](const std::vector<std::string>& v, char**& list, int& count) {
            count = 0;
            if (v.empty()) return;
            list = static_cast<char**>(checkedCalloc(v.size(), sizeof(char*)));
            for (const std::string& str : v) list[count++] = mallocString(str);
        };
        strings(aps, out.aps, out.noAPs);
        strings(accNameParameters, out.accNameParameters, out.noANPs);
        strings(properties, out.properties, out.noProps);

        out.start = mallocArray(start);
        out.noStart = static_cast<int>(start.size());
        out.cntAPs = mallocArray(cntAPs);
        out.noCntAPs = static_cast<int>(cntAPs.size());

        out.noAccSets = noAccSets;
        out.acc = acc.release();

        out.noAliases = 0;
        if (!aliases.empty()) {
            out.aliases = static_cast<Alias*>(checkedCalloc(aliases.size(), sizeof(Alias)));
            for (const AliasEntry& a : aliases) out.aliases[out.noAliases++].alias = mallocString(a.name);
        }

        out.noStates = 0;
        if (!states.empty()) {
            out.states = static_cast<State*>(checkedCalloc(states.size(), sizeof(State)));
            out.noStates = static_cast<int>(states.size());
            for (size_t i = 0; i < states.size(); i++) {
                out.states[i].id = static_cast<int>(i);
                out.states[i].accSig = mallocArray(states[i].accSig);
                out.states[i].noAccSig = static_cast<int>(states[i].accSig.size());
            }
        }
    } catch (...) {
        resetHoa(&out);
        throw;
    }

    // Then hand over the arena and the edges (cannot fail).
    for (size_t i = 0; i < aliases.size(); i++) out.aliases[i].labelExpr = aliases[i].expr;
    for (size_t i = 0; i < states.size(); i++) {
        State& dst = out.states[i];
        dst.name = states[i].name;
        dst.label = states[i].label;
        dst.transitions = states[i].transitions;
        dst.noTrans = states[i].noTrans;
    }
    out.arena = arena.release();
    stored = true;
    *data = out;
}

// Read the rest of the stream into a buffer followed by a NUL sentinel.
bool readAll(FILE* input, std::unique_ptr<char[]>& buffer, size_t& length)
{
    // For seekable streams, size the buffer so that one read reaches EOF.
    size_t capacity = size_t(1) << 16;
    const long here = ftell(input);
    if (here >= 0 && fseek(input, 0, SEEK_END) == 0) {
        const long last = ftell(input);
        if (fseek(input, here, SEEK_SET) != 0) return false;
        if (last >= here) capacity = static_cast<size_t>(last - here) + 2;
    }
    buffer.reset(new char[capacity]);
    length = 0;
    for (;;) {
        const size_t wanted = capacity - 1 - length;
        const size_t n = fread(buffer.get() + length, 1, wanted, input);
        length += n;
        if (n < wanted) break;
        std::unique_ptr<char[]> bigger(new char[capacity * 2]);
        memcpy(bigger.get(), buffer.get(), length);
        buffer = std::move(bigger);
        capacity *= 2;
    }
    buffer[length] = '\0';
    return ferror(input) == 0;
}

int parseText(const char* text, size_t length, HoaData* data)
{
    Lexer lex(text, text + length);
    try {
        bool aborted = false;
        for (;;) {
            try {
                lex.next();
                if (lex.kind == Tok::End) {
                    lex.error(lex.start, aborted ? "no automaton after \"--ABORT--\"" : "empty input");
                }
                AutomatonParser parser(lex);
                parser.parse();
                parser.store(data);
                return 0;
            } catch (const Aborted&) {
                aborted = true;
            }
        }
    } catch (const ParseError& e) {
        auto lc = lex.lineColumn(e.position);
        fprintf(stderr, "Parsing error: %s [line %ld, column %ld]\n", e.message.c_str(), lc.first, lc.second);
    } catch (const std::bad_alloc&) {
        fprintf(stderr, "Parsing error: out of memory\n");
    }
    return 1;
}

// ---------------------------------------------------------------------------
// Printing in the HOA format
// ---------------------------------------------------------------------------

void printString(FILE* out, const char* s)
{
    fputc('"', out);
    for (; *s != '\0'; s++) {
        if (*s == '"' || *s == '\\') fputc('\\', out);
        fputc(*s, out);
    }
    fputc('"', out);
}

// Print a label expression or acceptance condition such that parsing it
// gives the same tree: operands of & and | that are themselves | (or &
// as right operand of &) are parenthesised, as are binary operands of !.
void printExpr(FILE* out, const BTree* t)
{
    auto operand = [out](const BTree* c, bool parens) {
        if (parens) fputc('(', out);
        printExpr(out, c);
        if (parens) fputc(')', out);
    };
    switch (t->type) {
        case NT_BOOL: fputc(t->id ? 't' : 'f', out); break;
        case NT_AP:
        case NT_SET: fprintf(out, "%d", t->id); break;
        case NT_ALIAS: fprintf(out, "@%s", t->alias); break;
        case NT_NOT:
            fputc('!', out);
            operand(t->left, t->left->type == NT_AND || t->left->type == NT_OR);
            break;
        case NT_FIN:
        case NT_INF:
            fputs(t->type == NT_FIN ? "Fin(" : "Inf(", out);
            printExpr(out, t->left);
            fputc(')', out);
            break;
        case NT_AND:
            operand(t->left, t->left->type == NT_OR);
            fputs(" & ", out);
            operand(t->right, t->right->type == NT_OR || t->right->type == NT_AND);
            break;
        case NT_OR:
            operand(t->left, false);
            fputs(" | ", out);
            operand(t->right, t->right->type == NT_OR);
            break;
    }
}

void printSets(FILE* out, const int* sets, int count)
{
    if (sets == nullptr) return;
    fputs(" {", out);
    for (int i = 0; i < count; i++) fprintf(out, i == 0 ? "%d" : " %d", sets[i]);
    fputc('}', out);
}

}  // namespace

// ---------------------------------------------------------------------------
// Public interface
// ---------------------------------------------------------------------------

extern "C" {

void defaultsHoa(HoaData* data)
{
    data->aps = nullptr;
    data->noAPs = -1;
    data->accNameParameters = nullptr;
    data->noANPs = -1;
    data->properties = nullptr;
    data->noProps = -1;
    data->states = nullptr;
    data->noStates = -1;
    data->aliases = nullptr;
    data->noAliases = -1;
    data->start = nullptr;
    data->noStart = -1;
    data->cntAPs = nullptr;
    data->noCntAPs = -1;
    data->noAccSets = -1;
    data->acc = nullptr;
    data->version = nullptr;
    data->accNameID = nullptr;
    data->toolName = nullptr;
    data->toolVersion = nullptr;
    data->name = nullptr;
    data->arena = nullptr;
}

void resetHoa(HoaData* data)
{
    auto freeStrings = [](char** list, int count) {
        if (list == nullptr) return;
        for (int i = 0; i < count; i++) free(list[i]);
        free(list);
    };
    const bool inArena = data->arena != nullptr;

    free(data->version);
    free(data->accNameID);
    free(data->toolName);
    free(data->toolVersion);
    free(data->name);
    freeStrings(data->aps, data->noAPs);
    freeStrings(data->accNameParameters, data->noANPs);
    freeStrings(data->properties, data->noProps);
    free(data->start);
    free(data->cntAPs);
    freeTree(data->acc);
    if (data->aliases != nullptr) {
        for (int i = 0; i < data->noAliases; i++) {
            free(data->aliases[i].alias);
            if (!inArena) freeTree(data->aliases[i].labelExpr);
        }
        free(data->aliases);
    }
    if (data->states != nullptr) {
        for (int i = 0; i < data->noStates; i++) {
            State& s = data->states[i];
            free(s.accSig);
            for (int j = 0; j < s.noTrans; j++) free(s.transitions[j].accSig);
            if (inArena) continue;
            free(s.name);
            freeTree(s.label);
            for (int j = 0; j < s.noTrans; j++) {
                freeTree(s.transitions[j].label);
                free(s.transitions[j].successors);
            }
            free(s.transitions);
        }
        free(data->states);
    }
    delete static_cast<Arena*>(data->arena);
    defaultsHoa(data);
}

int parseHoaBuffer(const char* text, size_t length, HoaData* data)
{
    std::unique_ptr<char[]> buffer;
    try {
        buffer.reset(new char[length + 1]);
    } catch (const std::bad_alloc&) {
        fprintf(stderr, "Parsing error: out of memory\n");
        return 1;
    }
    memcpy(buffer.get(), text, length);
    buffer[length] = '\0';
    return parseText(buffer.get(), length, data);
}

int parseHoa(FILE* input, HoaData* data)
{
    std::unique_ptr<char[]> buffer;
    size_t length = 0;
    try {
        if (!readAll(input, buffer, length)) {
            perror("Parsing error: cannot read input");
            return 1;
        }
    } catch (const std::bad_alloc&) {
        fprintf(stderr, "Parsing error: out of memory\n");
        return 1;
    }
    return parseText(buffer.get(), length, data);
}

void printHoa(const HoaData* data)
{
    FILE* out = stdout;
    fprintf(out, "HOA: %s\n", data->version != nullptr ? data->version : "v1");
    if (data->name != nullptr) {
        fputs("name: ", out);
        printString(out, data->name);
        fputc('\n', out);
    }
    if (data->toolName != nullptr) {
        fputs("tool: ", out);
        printString(out, data->toolName);
        if (data->toolVersion != nullptr) {
            fputc(' ', out);
            printString(out, data->toolVersion);
        }
        fputc('\n', out);
    }
    if (data->noStates >= 0) fprintf(out, "States: %d\n", data->noStates);
    for (int i = 0; i < data->noStart; i++) fprintf(out, "Start: %d\n", data->start[i]);
    fprintf(out, "AP: %d", data->noAPs > 0 ? data->noAPs : 0);
    for (int i = 0; i < data->noAPs; i++) {
        fputc(' ', out);
        printString(out, data->aps[i]);
    }
    fputc('\n', out);
    if (data->noCntAPs >= 0) {
        fputs("controllable-AP:", out);
        for (int i = 0; i < data->noCntAPs; i++) fprintf(out, " %d", data->cntAPs[i]);
        fputc('\n', out);
    }
    for (int i = 0; i < data->noAliases; i++) {
        fprintf(out, "Alias: @%s ", data->aliases[i].alias);
        printExpr(out, data->aliases[i].labelExpr);
        fputc('\n', out);
    }
    if (data->accNameID != nullptr) {
        fprintf(out, "acc-name: %s", data->accNameID);
        for (int i = 0; i < data->noANPs; i++) fprintf(out, " %s", data->accNameParameters[i]);
        fputc('\n', out);
    }
    if (data->acc != nullptr) {
        fprintf(out, "Acceptance: %d ", data->noAccSets);
        printExpr(out, data->acc);
        fputc('\n', out);
    }
    if (data->noProps > 0) {
        fputs("properties:", out);
        for (int i = 0; i < data->noProps; i++) fprintf(out, " %s", data->properties[i]);
        fputc('\n', out);
    }
    fputs("--BODY--\n", out);
    for (int i = 0; data->states != nullptr && i < data->noStates; i++) {
        const State& s = data->states[i];
        fputs("State: ", out);
        if (s.label != nullptr) {
            fputc('[', out);
            printExpr(out, s.label);
            fputs("] ", out);
        }
        fprintf(out, "%d", s.id);
        if (s.name != nullptr) {
            fputc(' ', out);
            printString(out, s.name);
        }
        printSets(out, s.accSig, s.noAccSig);
        fputc('\n', out);
        for (int j = 0; j < s.noTrans; j++) {
            const Transition& t = s.transitions[j];
            if (t.label != nullptr) {
                fputc('[', out);
                printExpr(out, t.label);
                fputs("] ", out);
            }
            for (int k = 0; k < t.noSucc; k++) fprintf(out, k == 0 ? "%d" : "&%d", t.successors[k]);
            printSets(out, t.accSig, t.noAccSig);
            fputc('\n', out);
        }
    }
    fputs("--END--\n", out);
}


int isParityGFG(const HoaData* data, bool* isMaxParity, short* resGoodPriority)
{
    // (1) a parity acceptance condition in canonical form
    switch (parityName(data, isMaxParity, resGoodPriority)) {
        case 100:
            fprintf(stderr, "Not a parity automaton: acc-name is \"%s\", not \"parity\"\n",
                    data->accNameID != nullptr ? data->accNameID : "");
            return 100;
        case 101:
            fprintf(stderr, "The parity acc-name does not specify max or min\n");
            return 101;
        case 102:
            fprintf(stderr, "The parity acc-name does not specify even or odd\n");
            return 102;
        default:
            break;
    }
    const int first = *isMaxParity ? data->noAccSets - 1 : 0;
    if (!isParityCondition(data->acc, data->noAccSets, *isMaxParity, *resGoodPriority, first, false)) {
        fprintf(stderr, "The Acceptance: condition does not match its parity acc-name\n");
        return 103;
    }
    // (2) deterministic, complete and colored
    if (!hasString(data->properties, data->noProps, "deterministic")) {
        fprintf(stderr, "The automaton is not declared deterministic (properties: deterministic)\n");
        return 200;
    }
    if (!hasString(data->properties, data->noProps, "complete")) {
        fprintf(stderr, "The automaton is not declared complete (properties: complete)\n");
        return 201;
    }
    if (!hasString(data->properties, data->noProps, "colored")) {
        fprintf(stderr, "The automaton is not declared colored (properties: colored)\n");
        return 202;
    }
    // (3) a unique initial state
    if (data->noStart != 1) {
        fprintf(stderr, "The automaton must have exactly one initial state\n");
        return 300;
    }
    return 0;
}

}  // extern "C"
