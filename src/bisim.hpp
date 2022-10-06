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

TASK_DECL_2(sylvan::MTBDD, min_lts_strong, SymGame*, bool);
VOID_TASK_DECL_3(minimize, SymGame*, sylvan::MTBDD, bool);
VOID_TASK_DECL_2(print_partition, SymGame*, sylvan::MTBDD)
VOID_TASK_DECL_2(print_signature, SymGame*, sylvan::MTBDD)

size_t count_blocks();
