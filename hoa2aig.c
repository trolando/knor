/**************************************************************************
 * Copyright (c) 2019- Guillermo A. Perez
 * 
 * This file is part of HOATOOLS.
 * 
 * HOATOOLS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * HOATOOLS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with HOATOOLS. If not, see <http://www.gnu.org/licenses/>.
 * 
 * Guillermo A. Perez
 * University of Antwerp
 * guillermoalberto.perez@uantwerpen.be
 *************************************************************************/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "aiger/aiger.h"

#include "simplehoa.h"

/* A red-black tree to keep track of all and gates
 * we create for the automaton's transition function
 * Conventions:
 * (1) literals are positive or negative integers depending on whether they
 * are negated
 * (2) the key is composed of the left operand (literal) and the right operand
 * with the left operand being the "smaller" variable (not literal)
 * (3) variables are indexed from 2 (so 1 is "True" and -1 is "False")
 */
typedef enum {BLACK, RED} NodeColor;

typedef struct RBTree {
    struct RBTree* parent;
    struct RBTree* left;
    struct RBTree* right;
    NodeColor color;
    int opLeft;
    int opRight;
    int var;
} RBTree;

static char* colorStr(NodeColor c) {
    return c == RED ? "red" : "black";
}

static void recursivePrint(RBTree* n, int h) {
    // Essentially: a DFS with in-order printing
    if (n == NULL)
        return;
    recursivePrint(n->left, h + 1);
    fprintf(stderr, "%d: key=(%d,%d), var=%d (%s)\n", h, n->opLeft,
            n->opRight, n->var, colorStr(n->color));
    recursivePrint(n->right, h + 1);
}

static inline unsigned var2aiglit(int var) {
    if (var == -1)
        return 0;
    if (var == 1)
        return 1;
    int aiglit = 2 * (abs(var) - 1);
    if (var < 0)
        aiglit += 1;
    return aiglit;
}

static void dumpAiger(RBTree* n, aiger* aig) {
    // Essentially: a DFS with in-order dump
    if (n == NULL)
        return;
    dumpAiger(n->left, aig);
#ifndef NDEBUG
    fprintf(stderr, "Adding variable %d (%d)\n", var2aiglit(n->var), n->var);
    fprintf(stderr, "tis an and of %d (%d) with %d (%d)\n",
            var2aiglit(n->opLeft), n->opLeft,
            var2aiglit(n->opRight), n->opRight);
#endif
    aiger_add_and(aig,
                  var2aiglit(n->var),
                  var2aiglit(n->opLeft),
                  var2aiglit(n->opRight));
    dumpAiger(n->right, aig);
}

static void deleteTree(RBTree* n) {
    if (n == NULL)
        return;
    deleteTree(n->left);
    deleteTree(n->right);
    free(n);
}

static inline bool lessThan(RBTree* n, RBTree* m) {
    int nLeft = n->opLeft;
    int nRight = n->opRight;
    int mLeft = m->opLeft;
    int mRight = m->opRight;
    return (nLeft < mLeft) ||
           ((nLeft == mLeft)
            && (nRight < mRight));
}

static inline bool equalKeys(RBTree* n, RBTree* m) {
    return (n->opLeft == m->opLeft) && (n->opRight == m->opRight);
}

static RBTree* recursiveInsert(RBTree* root, RBTree* n) {
    // Essentially: go to the leaves and insert a red node
    if (root != NULL) {
        if (lessThan(n, root)) {
            if (root->left != NULL) {
                return recursiveInsert(root->left, n);
            } else {
                root->left = n;
            }
        } else if (equalKeys(n, root)) {
            return root;
        } else { // the root's key is strictly larger
            if (root->right != NULL) {
                return recursiveInsert(root->right, n);
            } else {
                root->right = n;
            }
        }
    }
    // insert n
    // NOTE: root may be NULL here
    n->parent = root;
    n->left = NULL;
    n->right = NULL;
    n->color = RED;
    return n;
}

static inline void attachToParent(RBTree* n, RBTree* p, RBTree* nnew) {
    if (p != NULL) {
        if (n == p->left) {
            p->left = nnew;
        } else if (n == p->right) {
            p->right = nnew;
        }
    }
    nnew->parent = p;
}

static void rotateLeft(RBTree* n) {
    RBTree* nnew = n->right;
    RBTree* p = n->parent;
    assert(nnew != NULL);  // Since the leaves of a red-black tree are empty,
                           // they cannot become internal nodes.
    n->right = nnew->left;
    nnew->left = n;
    n->parent = nnew;
    // handle other child/parent pointers
    if (n->right != NULL)
        n->right->parent = n;

    // attach to parent if n is not the root
    attachToParent(n, p, nnew);
}

static void rotateRight(RBTree* n) {
    RBTree* nnew = n->left;
    RBTree* p = n->parent;
    assert(nnew != NULL);  // Since the leaves of a red-black tree are empty,
                           // they cannot become internal nodes.
    n->left = nnew->right;
    nnew->right = n;
    n->parent = nnew;
    // handle other child/parent pointers
    if (n->left != NULL)
        n->left->parent = n;

    // attach to parent if n is not the root
    attachToParent(n, p, nnew);
}

static inline RBTree* sibling(RBTree* n) {
    if (n->parent == NULL) {
        return NULL;
    } else {
        return n == n->parent->left ? n->parent->right : n->parent->left;
    }
}

static inline RBTree* uncle(RBTree* n) {
    if (n->parent == NULL) {
        return NULL;
    } else {
        return sibling(n->parent);
    }
}

static void repairRBTree(RBTree* n) {
    if (n->parent == NULL) {
        // the root must be black
        n->color = BLACK;
    } else if (n->parent->color == BLACK) {
        // do nothing, all properties hold
    } else if (uncle(n) != NULL && uncle(n)->color == RED) {
        // we swap colors of parent & uncle with their parent
        n->parent->color = BLACK;
        uncle(n)->color = BLACK;
        n->parent->parent->color = RED;
        // and recurse on the grandparent to fix potential problems
        repairRBTree(n->parent->parent);
    } else {
        // we now want to rotate n with its grandparent, but for this to work
        // we need n to be either the left-left grandchild or the right-right
        // grandchild
        RBTree* p = n->parent;
        RBTree* g = p->parent;
        // we start by rotating n with its parent, if needed,
        // to guarantee this property
        if (n == p->right && p == g->left) {
            rotateLeft(p);
            n = n->left;
        } else if (n == p->left && p == g->right) {
            rotateRight(p);
            n = n->right;
        }
        // now we can rotate with the grandparent, and swap the colors of
        // parent and grandparent
        p = n->parent;
        g = p->parent;
        if (n == p->left) {
            rotateRight(g);
        } else {
            rotateLeft(g);
        }
        p->color = BLACK;
        g->color = RED;
    }
}

static RBTree* insertNode(RBTree* root, int op1, int op2, int freshVar) {
    RBTree* n;
    n = malloc(sizeof(RBTree));
    if (abs(op1) <= abs(op2)) {
        n->opLeft = op1;
        n->opRight = op2;
    } else {
        n->opLeft = op2;
        n->opRight = op1;
    }
    n->var = freshVar;
    n->parent = NULL;
    // the recursive insertion fails if a node with the
    // same key exists already
    RBTree* m = recursiveInsert(root, n);
    // in that case m is the node in the tree with the
    // same key but different variable
    if (m->var != n->var) {
        free(n);
        return m;
    }
    // repair the tree to recover red-black properties
    repairRBTree(n);
    return n;
}

static inline RBTree* findRoot(RBTree* n) {
    // find the new root
    RBTree* root = n;
    while (root->parent != NULL)
        root = root->parent;
    return root;
}

typedef struct {
    RBTree* root;
    int nextVar;
} AigTable;

static inline int and(AigTable* table, int op1, int op2) {
    assert(op1 != 0 && op2 != 0);
    RBTree* t = insertNode(table->root, op1, op2, table->nextVar);
    if (t->var == table->nextVar) {
        table->nextVar++;
        table->root = findRoot(t);
    }
    return t->var;
}

static inline int or(AigTable* table, int op1, int op2) {
    return -1 * and(table, -1 * op1, -1 * op2);
}

static int label2aig(AigTable* aig, BTree* label, AliasList* aliases) {
    assert(label != NULL);
    switch (label->type) {
        case NT_BOOL:
            return label->id ? 1 : -1;  // 0 becomes -1 like this
        case NT_AND:
            return and(aig, label2aig(aig, label->left, aliases),
                            label2aig(aig, label->right, aliases));
        case NT_OR:
            return or(aig, label2aig(aig, label->left, aliases),
                           label2aig(aig, label->right, aliases));
        case NT_NOT:
            return -1 * label2aig(aig, label->left, aliases);
        case NT_AP:
            return label->id + 2;  // FIXME: make this a global constant instead?
        case NT_ALIAS:
            for (AliasList* a = aliases; a != NULL; a = a->next) {
                if (strcmp(a->alias, label->alias) == 0)
                    return label2aig(aig, a->labelExpr, aliases);
            }
            break;
        default:
            assert(false);  // all cases should be covered above
    }
    return -2;
}

/* Read the EHOA file, encode the automaton in and-inverter
 * graphs, then use A. Biere's AIGER to dump the graph
 */
int main(int argc, char* argv[]) {
    HoaData* data = malloc(sizeof(HoaData));
    defaultsHoa(data);
    int ret = parseHoa(stdin, data);
    // 0 means everything was parsed correctly
    if (ret != 0)
        return ret;
    // A few semantic checks! TODO: use checkParityGFG function instead
    // (1) the automaton should be a parity one
    if (strcmp(data->accNameID, "parity") != 0) {
        fprintf(stderr, "Expected \"parity...\" automaton, found \"%s\" "
                        "as automaton type\n", data->accNameID);
        return 100;
    }
    bool foundOrd = false;
    bool maxPriority;
    bool foundRes = false;
    short winRes;
    for (StringList* param = data->accNameParameters; param != NULL;
            param = param->next) {
        if (strcmp(param->str, "max") == 0) {
            maxPriority = true;
            foundOrd = true;
        }
        if (strcmp(param->str, "min") == 0) {
            maxPriority = false;
            foundOrd = true;
        }
        if (strcmp(param->str, "even") == 0) {
            winRes = 0;
            foundRes = true;
        }
        if (strcmp(param->str, "odd") == 0) {
            winRes = 1;
            foundRes = true;
        }
    }
    if (!foundOrd) {
        fprintf(stderr, "Expected \"max\" or \"min\" in the acceptance name\n");
        return 101;
    }
    if (!foundRes) {
        fprintf(stderr, "Expected \"even\" or \"odd\" in the acceptance name\n");
        return 102;
    }
    // (2) the automaton should be deterministic, complete, colored
    bool det = false;
    bool complete = false;
    bool colored = false;
    for (StringList* prop = data->properties; prop != NULL; prop = prop->next) {
        if (strcmp(prop->str, "deterministic") == 0)
            det = true;
        if (strcmp(prop->str, "complete") == 0)
            complete = true;
        if (strcmp(prop->str, "colored") == 0)
            colored = true;
    }
    if (!det) {
        fprintf(stderr, "Expected a deterministic automaton, "
                        "did not find \"deterministic\" in the properties\n");
        return 200;
    }
    if (!complete) {
        fprintf(stderr, "Expected a complete automaton, "
                        "did not find \"complete\" in the properties\n");
        return 201;
    }
    if (!colored) {
        fprintf(stderr, "Expected one acceptance set per transition, "
                        "did not find \"colored\" in the properties\n");
        return 202;
    }
    // (3) the automaton should have a unique start state
    if (data->start == NULL || data->start->next != NULL) {
        fprintf(stderr, "Expected a unique start state\n");
        return 300;
    }

    // We now encode the transition relation into our
    // "sorta unique" AIG symbol table
    AigTable andGates;
    andGates.root = NULL;
    andGates.nextVar = 2;

    // We need to reserve a few variables though
    // (1) one per input + an extra input we will introduce
    // (2) one per latch needed for the states + 2 extra that we will need to
    // simulate eventual safety via safety + one per "good" acceptance set
    int noInputs = data->noAPs + 1;
    andGates.nextVar += noInputs;
    // we need floor(lg(#states)) + 1 just for the states,
    // where lg is the logarithm base 2; but C only has
    // natural logarithms
    int goodAccSets = (int) (ceil(data->noAccSets / 2.0));
    int noLatches = (int) (log(data->noStates) / log(2.0)) + 3 + goodAccSets;
    andGates.nextVar += noLatches;
#ifndef NDEBUG
    fprintf(stderr, "Reserved %d variables for %d inputs\n",
            noInputs, noInputs);
    fprintf(stderr, "Reserved %d variables for latches to encode %d states\n",
            noLatches, data->noStates);
    fprintf(stderr, "of them, %d are for good acceptance sets\n", goodAccSets);
#endif

    // We will traverse states, their transitions, bits set to 1 in the
    // successor via the transition, and add (i.e. logical or in place)
    // the transition
    // Step 1: set up state encodings
    int statecode[data->noStates];
    for (int i = 0; i < data->noStates; i++)
        statecode[i] = 1;
    // AIGER assumes 0 is the initial state, so we need to give the start
    // state the right latch valuation
    int start = data->start->i;
    int nextId = 1;
    int stateIds[data->noStates];
    for (int src = 0; src < data->noStates; src++) {
        // Step 1.1: create the encoding of the source state based on the
        // binary encoding of its number
        int mask = 1;
        int curId;
        if (src == start) {
            curId = 0;
        } else {
            curId = nextId;
            nextId++;
        }
        stateIds[src] = curId;
        for (int latch = 0; latch < noLatches - 2 - goodAccSets; latch++) {
            int latchlit = 2 + noInputs + latch;
            if ((curId & mask) != mask)
                latchlit *= -1;
            statecode[src] = and(&andGates, statecode[src], latchlit);
            mask = mask << 1;
        }
    }
    // Step 2: set up trivial predecessor and acceptance functions
    int predecessors[noLatches];
    for (int i = 0; i < noLatches; i++)
        predecessors[i] = -1;
    int acceptance[data->noAccSets];
    for (int i = 0; i < data->noAccSets; i++)
        acceptance[i] = -1;
    // Step 3: traverse the list of successors to update the latter
    for (StateList* src = data->states; src != NULL; src = src->next) {
        bool labelled = false;
        int labelCode;
        if (src->label != NULL) {
            labelCode = label2aig(&andGates, src->label, data->aliases);
            labelled = true;
        }

        for (TransList* trans = src->transitions; trans != NULL;
                                                  trans = trans->next) {
            if (!labelled) {
                if (trans->label == NULL) {
                    fprintf(stderr, "Cannot determine the label of a transition from state "
                                    "%d\n", src->id);
                    return 400;
                }
                labelCode = label2aig(&andGates, trans->label, data->aliases);
            }
            int noSuccessors = 0;
            for (IntList* tgt = trans->successors; tgt != NULL;
                                                   tgt = tgt->next) {
                int transCode = and(&andGates,
                                    labelCode,
                                    statecode[src->id]);
                int mask = 1;
                for (int latch = 0; latch < noLatches - 2 - goodAccSets; latch++) {
                    if ((stateIds[tgt->i] & mask) == mask) {
                        // Step 3.1: for each latch which should be set to 1
                        // for this successor, we update its predecessor
                        // relation
                        predecessors[latch] = or(&andGates, transCode,
                                                 predecessors[latch]);
                    }
                    mask = mask << 1;
                }
                noSuccessors++;
                // Step 3.2: we update the acceptance sets' functions based on
                // the transition and its acceptance sets
                IntList* acc = src->accSig;
                if (src->accSig == NULL)
                    acc = trans->accSig;
                // one of the two should be non-NULL
                // and there should be exactly one acceptance set!
                assert(acc != NULL && acc->next == NULL);
                acceptance[acc->i] = or(&andGates, acceptance[acc->i], transCode);
            }
            // since the automaton is supposedly deterministic, we should only
            // see one successor per transition
            assert(noSuccessors == 1);
        }
    }
    // Step 4: add logic regarding transitions of the sub-machine making sure
    // the new reset signal is used precisely once
    int resetInputVar = 2 + noInputs - 1;
    int resetLatch0 = noLatches - 2 - goodAccSets;
    int resetL0Var = 2 + noInputs + resetLatch0;
    int resetLatch1 = noLatches - 1 - goodAccSets;
    int resetL1Var = 2 + noInputs + resetLatch1;
    // we move to accepting (i.e. ~resetLatch1 & resetLatch0) if
    // we are at the initial state and read a reset input
    int initialState = and(&andGates, -1 * resetL0Var, -1 * resetL1Var);
    predecessors[resetLatch0] = or(&andGates, predecessors[resetLatch0],
                                   and(&andGates, resetInputVar, initialState));
    // or if we are already accepting and read no reset input
    int acceptState = and(&andGates, -1 * resetL1Var, resetL0Var);
    predecessors[resetLatch0] = or(&andGates, predecessors[resetLatch0],
                                   and(&andGates, -1 * resetInputVar, acceptState));
    // we move to losing (i.e. resetLatch1 & ~resetLatch0) if
    // we are at the accepting state and read a reset input
    predecessors[resetLatch1] = or(&andGates, predecessors[resetLatch1],
                                   and(&andGates, resetInputVar, acceptState));
    // or if we are already losing and read anything
    int loseState = and(&andGates, resetL1Var, -1 * resetL0Var);
    predecessors[resetLatch1] = or(&andGates, predecessors[resetLatch1],
                                   loseState);
    // Step 5: we prepare the logic for the justice conditions based on the type
    // of parity condition
    // NOTE: winRes = 1 for odd parity conditions, winRes = 0 otherwise
    int safeResetLatch = noLatches - goodAccSets;
    for (int p = 1 - winRes; p < data->noAccSets; p += 2) {
        int trumpP = -1;
        // NOTE: maxPriority = true iff we have a max priority condition
        int q = maxPriority ? p + 1 : p - 1;
        while (q >= 0 && q < data->noAccSets) {
            trumpP = or(&andGates, trumpP, acceptance[q]);
            q += maxPriority ? 2 : -2;
        }
        int safeResetVar = 2 + noInputs + safeResetLatch;
        predecessors[safeResetLatch] = or(&andGates, safeResetVar,
                                          and(&andGates,
                                              acceptState,
                                              trumpP));
        safeResetLatch++;
    }

#ifndef NDEBUG
    recursivePrint(andGates.root, 0);
#endif

    // Step 6: Create and print the constructed AIG
    aiger* aig = aiger_init();
    // add inputs
    int lit = 2;
    for (StringList* aps = data->aps; aps != NULL; aps = aps->next) {
        aiger_add_input(aig, lit, aps->str);
        lit += 2;
    }
    // add the extra input
    aiger_add_input(aig, lit, "reset_safety");
    lit += 2;
    // add latches
    for (int latch = 0; latch < noLatches - 2 - goodAccSets; latch++) {
        aiger_add_latch(aig, lit, var2aiglit(predecessors[latch]), ""); 
        lit += 2;
    }
    aiger_add_latch(aig, lit,
                    var2aiglit(predecessors[noLatches - 2 - goodAccSets]),
                    "reset_once_latch0");
    lit += 2;
    aiger_add_latch(aig, lit,
                    var2aiglit(predecessors[noLatches - 1 - goodAccSets]),
                    "reset_once_latch1");
    lit += 2;
    char latchName[50];
    for (int latch = noLatches - goodAccSets; latch < noLatches; latch++) {
        int localNo = latch - (noLatches - goodAccSets);
        sprintf(latchName, "safe_with_1_reset%d", localNo);
        aiger_add_latch(aig, lit, var2aiglit(predecessors[latch]), latchName);
        lit += 2;
    }
    
    // add and-gates
#ifndef NDEBUG
    fprintf(stderr, "Dumping AND-gates into aiger structure\n");
#endif
    dumpAiger(andGates.root, aig);
    // add fairness constraint for the extra input
    aiger_add_fairness(aig, var2aiglit(acceptState), "reset_safety_once");
    // add justice constraints
    safeResetLatch = noLatches - goodAccSets;
    char justiceName[50];
    for (int p = 1 - winRes; p < data->noAccSets; p += 2) {
        int safeResetVar = 2 + noInputs + safeResetLatch;
        unsigned justice[2] = {var2aiglit(-1 * safeResetVar),
                               var2aiglit(acceptance[p])};
        sprintf(justiceName, "priority_%d_untrumped", p);
        aiger_add_justice(aig, 2, justice, justiceName);
        safeResetLatch++;
    }

#ifndef NDEBUG
    fprintf(stderr, "AIG structure created, now checking it!\n");
    const char* msg = aiger_check(aig);
    if (msg) {
        fprintf(stderr, msg);
        return 500;
    }
#endif

    // and dump the aig
    aiger_write_to_file(aig, aiger_ascii_mode, stdout);

    // Free dynamic memory
    aiger_reset(aig);
    deleteTree(andGates.root);
    deleteHoa(data);
    return EXIT_SUCCESS;
}
