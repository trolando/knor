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

#include "simplehoa.h"

/* A red-black tree to keep track of all and gates
 * we create for the automaton's transition function
 * Conventions:
 * (1) literals are positive or negative integers
 *     depending on whether they are negated
 * (2) the key is composed of the left operand (literal)
 *     and the right operand with the left operand being
 *     the "smaller" variable (not literal)
 * (3) variables are indexed from 2 (so 1 is "True" and
 *     -1 is "False")
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
    printf("%d: key=(%d,%d), var=%d (%s); ", h, n->opLeft, n->opRight,
           n->var, colorStr(n->color));
    recursivePrint(n->right, h + 1);
}
#endif

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

/* Read the EHOA file, encode the automaton in and-inverter
 * graphs, then use A. Biere's AIGER to write the graph
 */
int main(int argc, char* argv[]) {
    HoaData* data = malloc(sizeof(HoaData));
    defaultsHoa(data);
    int ret = parseHoa(stdin, data);
    // 0 means everything was parsed correctly
    if (ret != 0)
        return ret;

    // We now encode the transition relation into our
    // "unique" AIG symbol table
    AigTable andGates;
    andGates.root = NULL;
    andGates.nextVar = 2;
    // I need to reserve a few variables though
    // (1) one per input + an extra input I will introduce
    // (2) one per latch needed for the states + 2
    //     extra that I will need to simulate FG
    andGates.nextVar += data->noAPs + 1;
    // we need floor(lg(#states)) + 1 just for the states,
    // where lg is the logarithm base 2; but C only has
    // natural logarithms
    int noLatches = (int) (log(data->noStates) / log(2.0)) + 1;
    andGates.nextVar += noLatches + 2;
#ifndef NDEBUG
    printf("Reserved %d variables for inputs\n",
           data->noAPs + 1);
    printf("Reserved %d variables for latches\n",
           noLatches + 2);
#endif
    int v = and(&andGates, 4, -2);
    v = and(&andGates, 5, 4);
    v = and(&andGates, 4, v * -1);
    //assert(v < andGates.nextVar);

#ifndef NDEBUG
    recursivePrint(andGates.root, 0);
    printf("\n");
#endif

    deleteTree(andGates.root);
    deleteHoa(data);
    return EXIT_SUCCESS;
}
