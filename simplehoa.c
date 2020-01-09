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

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>

#include "simplehoa.h"

IntList* newIntNode(int val) {
    IntList* list = malloc(sizeof(IntList));
    list->i = val;
    list->next = NULL;
    return list;
}

IntList* appendIntNode(IntList* node, int val) {
    IntList* newHead = malloc(sizeof(IntList));
    newHead->i = val;
    newHead->next = node;
    return newHead;
}

StringList* appendStrNode(StringList* node, char* str) {
    StringList* newHead = malloc(sizeof(StringList));
    newHead->str = malloc(sizeof(char) * strlen(str) + 1);
    strcpy(newHead->str, str);
    newHead->next = node;
    return newHead;
}

StringList* concatStrLists(StringList* list1, StringList* list2) {
    if (list2 == NULL)
        return list1;
    if (list1 == NULL)
        return list2;

    StringList* cur = list2;
    while (cur->next != NULL)
        cur = cur->next;
    cur->next = list1;
    return list2;
}

BTree* boolBTree(bool b) {
    BTree* created = malloc(sizeof(BTree));
    created->left = NULL;
    created->right = NULL;
    created->type = NT_BOOL;
    created->set = b ? 1 : 0;
    return created;
}

BTree* andBTree(BTree* u, BTree* v) {
    BTree* created = malloc(sizeof(BTree));
    created->left = u;
    created->right = v;
    created->type = NT_AND;
    created->set = -1;
    return created;
}

BTree* orBTree(BTree* u, BTree* v) {
    BTree* created = malloc(sizeof(BTree));
    created->left = u;
    created->right = v;
    created->type = NT_OR;
    created->set = -1;
    return created;
}

BTree* idBTree(NoteType type, int set, bool negated) {
    BTree* set = malloc(sizeof(BTree));
    set->left = NULL;
    set->right = NULL;
    set->type = NT_SET;
    set->set = set;

    if (negated) {
        BTree* original = set;
        set = malloc(sizeof(BTree));
        set->left = original;
        set->right = NULL;
        set->type = NT_NOT;
        set->set -1;
    }

    BTree* created = malloc(sizeof(BTree));
    created->left = set;
    created->right = NULL;
    created->type = type;
    created->set = -1;
    return created;
}

void defaultsHoa(HoaData* data) {
    data->noStates = -1;  // to say we have not gotten
    data->noAPs = -1;     // a number for these parameters
    data->start = NULL;
    data->version = NULL;
    data->aps = NULL;
    // data->noAccSets  // these need no default as they will
    // data->acc        // always be set by the parser
    data->accNameID = NULL;
    data->accNameParameters = NULL;
    data->toolName = NULL;
    data->toolVersion = NULL;
    data->name = NULL;
    data->properties = NULL;
    data->states = NULL;
}

static void deleteStrList(StringList* list) {
    StringList* cur = list;
    while (cur != NULL) {
        StringList* next = cur->next;
        if (cur->str != NULL)
            free(cur->str);
        free(cur);
        cur = next;
    }
}

static void deleteIntList(IntList* list) {
    IntList* cur = list;
    while (cur != NULL) {
        IntList* next = cur->next;
        free(cur);
        cur = next;
    }
}

void deleteHoa(HoaData* data) {
    // Strings
    if (data->version != NULL) free(data->version);
    if (data->accNameID != NULL) free(data->accNameID);
    if (data->toolName != NULL) free(data->toolName);
    if (data->toolVersion != NULL) free(data->toolVersion);
    if (data->name != NULL) free(data->name);
    // String lists
    deleteStrList(data->aps);
    deleteStrList(data->accNameParameters);
    deleteStrList(data->properties);
    // Int lists
    deleteIntList(data->start);
    // BTrees
    //BTree* acc;
    // State lists
    //StateList* staes;
}

void printHoa(const HoaData* data) {
    printf("HOA format version: %s\n", data->version);
    if (data->name != NULL)
        printf("File name: %s\n", data->name);
    printf("No. of states: %d\n", data->noStates);

    printf("Start states: ");
    for (IntList* it = data->start; it != NULL; it = it->next)
        printf("%d, ", it->i);
    printf("\n");

    printf("No. of atomic propositions: %d\n", data->noAPs);

    printf("Atomic propositions: ");
    for (StringList* it = data->aps; it != NULL; it = it->next)
        printf("%s, ", it->str);
    printf("\n");

    printf("No. of acceptance sets: %d\n", data->noAccSets);
    // TODO: print data->acc
    if (data->accNameID != NULL)
        printf("Acceptance name: %s\n", data->accNameID);

    printf("Acceptance parameters: ");
    for (StringList* it = data->accNameParameters; it != NULL; it = it->next)
        printf("%s, ", it->str);
    printf("\n");

    if (data->toolName != NULL) {
        assert(data->toolVersion != NULL);
        printf("Tool name: %s version %s\n",
               data->toolName,
               data->toolVersion);
    }

    printf("Properties: ");
    for (StringList* it = data->properties; it != NULL; it = it->next)
        printf("%s, ", it->str);
    printf("\n");
    // TODO: print data->states
}
