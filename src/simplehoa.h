/**
 * Copyright Tom van Dijk
 */

/*
 * In-memory representation of an automaton in the extended HOA format (eHOA):
 * the Hanoi Omega-Automata format v1 (https://adl.github.io/hoaf/) plus the
 * header item "controllable-AP:", which lists the indices of the atomic
 * propositions that are controlled by the system.
 *
 * This header is plain C so it can be included from C and C++ alike.
 */

#ifndef KNOR_SIMPLEHOA_H
#define KNOR_SIMPLEHOA_H

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Node types of label expressions and acceptance conditions. */
typedef enum {
    NT_BOOL,    /* constant: id is 1 (t) or 0 (f) */
    NT_AND,     /* left & right */
    NT_OR,      /* left | right */
    NT_FIN,     /* Fin(left) */
    NT_INF,     /* Inf(left) */
    NT_NOT,     /* !left */
    NT_SET,     /* acceptance set with number id */
    NT_AP,      /* atomic proposition with number id */
    NT_ALIAS    /* alias with name alias (without the '@') */
} NodeType;

/*
 * Expression tree node. Unused children are NULL, alias is NULL except for
 * NT_ALIAS, and id is -1 for the operators NT_AND, NT_OR, NT_NOT, NT_FIN and
 * NT_INF. Binary operators are left-associative: "a & b & c" is
 * AND(AND(a, b), c), and & binds stronger than |.
 */
typedef struct BTree {
    struct BTree* left;
    struct BTree* right;
    char* alias;
    int id;
    NodeType type;
} BTree;

typedef struct Transition {
    BTree* label;       /* NULL iff the source state carries the label */
    int* successors;    /* conjunction of successors (one unless alternating) */
    int noSucc;
    int* accSig;        /* acceptance sets, or NULL */
    int noAccSig;
} Transition;

typedef struct State {
    int id;             /* always equal to the index in HoaData.states */
    char* name;         /* or NULL */
    BTree* label;       /* or NULL */
    int* accSig;        /* acceptance sets, or NULL */
    int noAccSig;
    Transition* transitions;
    int noTrans;
} State;

typedef struct Alias {
    char* alias;        /* name without the '@' */
    BTree* labelExpr;
} Alias;

/*
 * A parsed automaton. Counts are -1 after defaultsHoa(); after a successful
 * parse, every list is NULL when its count is 0, and lists appear in file
 * order. Notable conventions:
 *  - states is indexed by state number and has noStates entries, where
 *    noStates is the "States:" value, or one more than the largest state
 *    number used if that header is absent. States without a "State:" entry
 *    in the body have no transitions.
 *  - all "Start:" items are concatenated into start; a conjunction of
 *    initial states (alternating automata) is flattened into the list.
 *  - acc-name parameters are stored as strings; t and f become "True" and
 *    "False", integers are written in decimal.
 *  - edges without a label in a state without a label (implicit labels)
 *    receive the label that the format assigns to them: edge i is labelled
 *    with the cube over all propositions that corresponds to the bits of i.
 *
 * Ownership. Everything is allocated with malloc and released by resetHoa,
 * so code that builds or edits an HoaData by hand may free and replace
 * fields with free/malloc, with one exception: when arena != NULL (as set
 * by the parser), the label trees of states, transitions and aliases, the
 * transitions arrays, the successors arrays and the state names live in
 * that arena. These are read-only and must not be freed individually; label
 * trees are shared (identical subexpressions are one node). The acceptance
 * condition acc, all accSig arrays, the states and aliases arrays and all
 * other strings and lists are always individually malloc'd.
 */
typedef struct HoaData {
    char** aps;
    int noAPs;
    char** accNameParameters;
    int noANPs;
    char** properties;
    int noProps;
    State* states;
    int noStates;
    Alias* aliases;
    int noAliases;
    int* start;
    int noStart;
    int* cntAPs;
    int noCntAPs;
    int noAccSets;
    BTree* acc;
    char* version;
    char* accNameID;
    char* toolName;
    char* toolVersion;
    char* name;
    void* arena;
} HoaData;

/*
 * Parse one automaton from the stream (a regular file, pipe, stdin or
 * fmemopen buffer are all fine; the stream is read until EOF). The data must
 * be in the state left by defaultsHoa(). Returns 0 on success; otherwise an
 * error message is printed to stderr, a nonzero value is returned and data
 * is left as it was. The input may be a stream of automata: one ended by
 * --ABORT-- is skipped, and the first complete one is returned (the ones
 * after it are not parsed).
 */
int parseHoa(FILE* input, HoaData* data);

/* As parseHoa, but parse the given buffer of length bytes. */
int parseHoaBuffer(const char* text, size_t length, HoaData* data);

/* Initialise all fields to "absent" (counts -1, pointers NULL). */
void defaultsHoa(HoaData* data);

/* Free everything owned by data and reinitialise it with defaultsHoa(). */
void resetHoa(HoaData* data);

/* Print the automaton to stdout in the (extended) HOA format. */
void printHoa(const HoaData* data);

/*
 * Rewrite the acceptance condition into an equivalent "parity max even"
 * condition with exactly one acceptance set per edge. Without memory, this
 * works for every parity condition in disguise, i.e., one that a priority per
 * edge expresses on the same automaton (its Zielonka tree over the
 * combinations of acceptance sets on edges is a chain): t and f, Buchi and
 * co-Buchi, Rabin and Streett with one pair, parity conditions in any order
 * or notation, edges in zero or several acceptance sets, and state-based
 * acceptance. With addMemory, any acceptance condition (generalized Buchi,
 * Rabin, Streett, Muller, ...) is converted, by replacing the automaton with
 * its product with the Zielonka tree automaton of the condition (this needs
 * an automaton read by the parser, without alternation).
 * An automaton that already has a canonical parity condition (as checked by
 * isParityGFG) with one acceptance set per edge or per state is left as it is,
 * except that the property "colored" is added if missing. Returns 0 on
 * success; otherwise prints the reason to stderr and returns 104 if the
 * condition needs memory but addMemory is false, or 100 otherwise.
 */
int paritizeHoa(HoaData* data, bool addMemory);

/*
 * Check that the automaton is a complete deterministic colored parity
 * automaton with a unique initial state, and determine its parity
 * condition. Returns 0 if so; otherwise prints the reason to stderr and
 * returns an error code: 100-103 (acceptance), 200-202 (properties) or
 * 300 (initial state).
 */
int isParityGFG(const HoaData* data, bool* isMaxParity, short* resGoodPriority);

#ifdef __cplusplus
}
#endif

#endif
