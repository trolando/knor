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

typedef enum {
    NT_BOOL, NT_AND, NT_OR, NT_FIN, NT_INF, NT_NOT
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
    IntList* start;
    char* version;
    int noAccSets;
    BTree* acc;
    char* accNameID;
    StringList* accNameParameters;
    char* toolName;
    char* toolVersion;
    char* name;
    StringList* properties;
    StateList* states;
} HoaData;

// This function returns 0 if and only if parsing was successful
int parseHoa(FILE*, HoaData*);

#endif
