/**
 * Copyright Tom van Dijk
 */

#include <sylvan.h>
#include <oink/oink.hpp>
#include <symgame.hpp>

extern "C" {
    #include "simplehoa.h"
}

#pragma once

// Bisimulation minimisation

TASK_DECL_1(sylvan::MTBDD, min_lts_strong, SymGame*);
VOID_TASK_DECL_3(minimize, SymGame*, sylvan::MTBDD, bool);
VOID_TASK_DECL_2(print_partition, SymGame*, sylvan::MTBDD)

