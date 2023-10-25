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
 * Encode priority i.e. all states via priority <priority>
 */
sylvan::MTBDD encode_prio(int priority, int priobits);

/**
 * Encode a state as a BDD, using statebits 0..<statebits>, offsetted by <offset>+<priobits>
 * High-significant bits come before low-significant bits in the BDD
 */
sylvan::MTBDD encode_state(uint32_t state, const int statebits, const int s_first_var);

/**
* Encode a priostate as a BDD, using priobits > statebits
* High-significant bits come before low-significant bits in the BDD
*/
sylvan::MTBDD encode_priostate(uint32_t state, uint32_t priority, const int statebits, const int priobits, const int offset);

/**
 * Convert a transition label (Btree) to a BDD encoding the label
 * a label is essentially a boolean combination of atomic propositions and aliases
 */
TASK_DECL_3(sylvan::MTBDD, evalLabel, BTree* /*label*/, HoaData* /*data*/, uint32_t* /*variables*/);
