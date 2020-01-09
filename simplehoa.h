/**************************************************************************
 * Copyright (c) 2019- Guillermo A. Perez
 * 
 * This file is part of the HOA2AIG tool.
 * 
 * HOA2AIG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * HOA2AIG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with HOA2AIG. If not, see <http://www.gnu.org/licenses/>.
 * 
 * Guillermo A. Perez
 * University of Antwerp
 * guillermoalberto.perez@uantwerpen.be
 *************************************************************************/

#ifndef _SIMPLEHOA_H
#define _SIMPLEHOA_H

#include <stdio.h>
#include <stdbool.h>

typedef enum {
    NT_BOOL, NT_AND, NT_OR, NT_FIN, NT_INF, NT_NOT, NT_SET
} NodeType;

typedef struct BTree BTree;
struct BTree {
    BTree* left;
    BTree* right;
    int set;
    NodeType type;
};

typedef struct StringList StringList;
struct StringList {
    char* str;
    StringList* next;
};

typedef struct IntList IntList;
struct IntList {
    int i;
    IntList* next;
};

typedef struct TransList TransList;
struct TransList {
    BTree* label;
    IntList* successors;
    IntList* accSig;
};

typedef struct StateList StateList;
struct StateList {
    int id;
    char* name;
    BTree* label;
    IntList* accSig;
    TransList* transitions;
};

typedef struct HoaData {
    int noStates;
    StringList* aps;
    StringList* accNameParameters;
    StringList* properties;
    StateList* states;
    IntList* start;
    int noAccSets;
    int noAPs;
    BTree* acc;
    char* version;
    char* accNameID;
    char* toolName;
    char* toolVersion;
    char* name;
} HoaData;

// This function returns 0 if and only if parsing was successful
int parseHoa(FILE*, HoaData*);
void defaultsHoa(HoaData*);
void deleteHoa(HoaData*);
IntList* newIntNode(int);
IntList* appendIntNode(IntList*, int);
StringList* appendStrNode(StringList*, char*);
StringList* concatStrLists(StringList*, StringList*);
BTree* boolBTree(bool);
BTree* andBTree(BTree*, BTree*);
BTree* orBTree(BTree*, BTree*);
BTree* idBTree(int, bool);

// For debugging purposes, this prints all data in human-readable form
void printHoa(const HoaData*);

#endif
