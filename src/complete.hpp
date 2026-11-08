/**
 * Copyright Tom van Dijk
 */

#ifndef KNOR_COMPLETE_HPP
#define KNOR_COMPLETE_HPP

#include "simplehoa.h"

/*
 * Check that the automaton is deterministic and complete, as far as it does
 * not declare these properties, and declare them. An automaton read by the
 * parser that is not complete is completed: the valuations without an edge
 * lead to a new rejecting sink state (see complete.cpp for its acceptance).
 * The labels are compared as BDDs, with proposition i as variable i, so
 * Sylvan must be initialised. Returns 0 on success; otherwise prints the
 * reason to stderr and returns 200 if the automaton is not deterministic, or
 * 100 otherwise.
 */
int completeHoa(HoaData* data);

#endif
