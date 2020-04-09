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

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "aiger.h"

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

typedef struct RBTree RBTree;
struct RBTree {
    RBTree* parent;
    RBTree* left;
    RBTree* right;
    NodeColor color;
    int opLeft;
    int opRight;
    int var;
};

#ifndef NDEBUG
static char* colorStr(NodeColor c) {
    return c == RED ? "red" : "black";
}

static void recursivePrint(RBTree* n, int h) {
    // Essentially: a DFS with in-order printing
    if (n == NULL)
        return;
    recursivePrint(n->left, h + 1);
    printf("%d: key=(%d,%d), var=%d (%s)\n", h, n->opLeft, n->opRight,
           n->var, colorStr(n->color));
    recursivePrint(n->right, h + 1);
}
#endif

static void dumpAiger(RBTree* n, aiger* aig) {
    // Essentially: a DFS with in-order dump
    if (n == NULL)
        return;
    dumpAiger(n->left, aig);
    unsigned var = 2 * (n->var - 1);
    unsigned op1 = 2 * (abs(n->opLeft) - 1);
    if (n->opLeft < 0)
        op1 += 1;
    unsigned op2 = 2 * (abs(n->opRight) - 1);
    if (n->opRight < 0)
        op2 += 1;
    aiger_add_and(aig, var, op1, op2);
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
    // A few semantic checks!
    // (1) the automaton should be a parity one
    if (strncmp(data->accNameID, "parity", 6) != 0) {
        fprintf(stderr, "Expected \"parity...\" automaton, found \"%s\" "
                        "as automaton type\n", data->accNameID);
        return 100;
    }
    // (2) the automaton should be deterministic
    bool det = false;
    for (StringList* prop = data->properties; prop != NULL; prop = prop->next) {
        if (strncmp(prop->str, "deterministic", 13) == 0)
            det = true;
    }
    if (!det) {
        fprintf(stderr, "Expected a deterministic automaton, "
                        "did not find \"deterministic\" in the properties\n");
        return 200;
    }

    // We now encode the transition relation into our
    // "sorta unique" AIG symbol table
    AigTable andGates;
    andGates.root = NULL;
    andGates.nextVar = 2;

    // We need to reserve a few variables though
    // (1) one per input + an extra input we will introduce
    // (2) one per latch needed for the states + 2 extra that we will need to
    // simulate eventual safety via safety
    int noInputs = data->noAPs + 1;
    andGates.nextVar += noInputs;
    // we need floor(lg(#states)) + 1 just for the states,
    // where lg is the logarithm base 2; but C only has
    // natural logarithms
    int noLatches = (int) (log(data->noStates) / log(2.0)) + 3;
    andGates.nextVar += noLatches;
#ifndef NDEBUG
    printf("Reserved %d variables for inputs\n",
           noInputs);
    printf("Reserved %d variables for latches to encode %d states\n",
           noLatches, data->noStates);
#endif

    // We will traverse states, their transitions, bits set to 1 in the
    // successor via the transition, and add (i.e. logical or in place)
    // the transition
    // Step 1: set up state encodings
    int statecode[data->noStates];
    for (int i = 0; i < data->noStates; i++)
        statecode[i] = 1;
    for (int src = 0; src < data->noStates; src++) {
        // Step 1.1: create the encoding of the source state based on the
        // binary encoding of its number
        int mask = 1;
        for (int latch = 0; latch < noLatches - 2; latch++) {
            int latchlit = noInputs + latch;
            if ((src & mask) != mask)
                latchlit *= -1;
            statecode[src] = and(&andGates, statecode[src], latchlit);
            mask = mask << 1;
        }
    }
    // Step 2: set up trivial predecessor relations
    int predecessors[noLatches];
    memset(predecessors, 0, noLatches * sizeof(int));
    // Step 3: traverse the list of successors to update the latter
    for (StateList* src = data->states; src != NULL; src = src->next) {
        for (TransList* trans = src->transitions; trans != NULL;
                                                  trans = trans->next) {
            int noSuccessors = 0;
            for (IntList* tgt = trans->successors; tgt != NULL;
                                                   tgt = tgt->next) {
                int mask = 1;
                for (int latch = 0; latch < noLatches - 2; latch++) {
                    if (tgt->i & mask == 1) {
                        // Step 3.1: for each latch which should be set to 1
                        // for this successor, we update its predecessor
                        // relation
                        // FIXME: the transition label on the inputs is
                        // missing
                        predecessors[latch] = or(&andGates,
                                                 predecessors[latch],
                                                 statecode[src->id]);
                    }
                    mask << 1;
                }
                noSuccessors++;
            }
            // since the automaton is supposedly deterministic, we should only
            // see one successor per transition
            assert(noSuccessors == 1);
        }
    }

#ifndef NDEBUG
    recursivePrint(andGates.root, 0);
#endif
    // Create and print the constructed AIG
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
    char latchName[50];
    for (int latch = 0; latch < noLatches; latch++) {
        sprintf("latch_%d", latchName, latch);
        aiger_add_latch(aig, lit, 2 * (predecessors[latch] - 1), latchName); 
        lit += 2;
    }
    // add and-gates
    dumpAiger(andGates.root, aig);
#ifndef NDEBUG
    aiger_check(aig);
#endif
    // and dump the aig
    aiger_write_to_file(aig, aiger_ascii_mode, stdout);

    // Free dynamic memory
    aiger_reset(aig);
    deleteTree(andGates.root);
    deleteHoa(data);
    return EXIT_SUCCESS;
}
