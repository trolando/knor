/**
 * Copyright Tom van Dijk
 */

#include <sylvan.h>
#include <oink/oink.hpp>

extern "C" {
#include "simplehoa.h"
}

#pragma once

/**
 * This helper function ensures that the priority p is adjusted to
 * ensure we have a "max, even" parity game, as this is what Oink expects.
 * @param controllerIsOdd if the controller is the odd player
 * @param noPriorities how many priorities are in the game
 * @param maxPriority if the game is a max game
 */
int adjustPriority(int p, bool maxPriority, bool controllerIsOdd, int noPriorities);

/**
 * Convert a transition label (Btree) to a BDD encoding the label
 * a label is essentially a boolean combination of atomic propositions and aliases
 */
TASK_DECL_3(sylvan::MTBDD, evalLabel, BTree* /*label*/, HoaData* /*data*/, uint32_t* /*variables*/);
