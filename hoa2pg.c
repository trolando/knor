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
#include <string.h>

#include "cudd/cudd.h"

#include "simplehoa.h"


/* Read the EHOA file, construct a graph-based game and
 * dump it as a PGSolver game
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
    
    printf("Everything loaded successfully!\n");

    // TODO: do something here

    // Free dynamic memory
    deleteHoa(data);
    return EXIT_SUCCESS;
}
