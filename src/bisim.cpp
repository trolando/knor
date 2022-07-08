/**
 * Copyright Tom van Dijk
 */

#include <vector>
#include <algorithm>
#include <iostream>
#include <sylvan.h>
#include <sylvan_int.h>
#include <sys/mman.h>

#include <knor.hpp>
#include <bisim.hpp>

using namespace sylvan;

int block_length;                 // number of block variables
MTBDD block_variables;            // BDD of block variables

uint64_t CACHE_ENCODE_BLOCK = 0;
uint64_t CACHE_DECODE_BLOCK = 0;
uint64_t CACHE_REFINE = 0;
uint64_t CACHE_APPLY_PARTITION = 0;

TASK_1(BDD, encode_block, uint64_t, b)
{
    BDD result;
    if (cache_get3(CACHE_ENCODE_BLOCK, 0, b, 0, &result)) return result;

    // for now, assume max 64 bits for a block....
    uint8_t bl[block_length];
    for (int i=0; i<block_length; i++) {
        bl[i] = b & 1 ? 1 : 0;
        b>>=1;
    }

    result = sylvan_cube(block_variables, bl);
    cache_put3(CACHE_ENCODE_BLOCK, 0, b, 0, result);
    return result;
}

TASK_1(uint64_t, decode_block, BDD, block)
{
    uint64_t result = 0;
    if (cache_get3(CACHE_DECODE_BLOCK, block, 0, 0, &result)) return result;

    uint64_t mask = 1;
    while (block != sylvan_true) {
        BDD b_low = sylvan_low(block);
        if (b_low == sylvan_false) {
            result |= mask;
            block = sylvan_high(block);
        } else {
            block = b_low;
        }
        mask <<= 1;
    }

    cache_put3(CACHE_DECODE_BLOCK, block, 0, 0, result);
    return result;
}

/**
 * Signatures records the signature of each block.
 * The default and maximum size of this array is "signatures_size", 1LL<<25
 * In practice, we will run out of memory and time way before this.
 */

#define SL_DEPTH 5
typedef struct
{
    BDD sig;
    uint32_t prev;
    uint32_t next[SL_DEPTH];
} signature_elem;

static signature_elem *signatures = NULL;
static size_t signatures_size = 0;
static uint32_t next_block = 1;
static size_t refine_iteration = 0;

void
prepare_refine()
{
    if (signatures != NULL) {
        munmap(signatures, sizeof(signature_elem)*signatures_size);
    }
    signatures = (signature_elem*)mmap(0, sizeof(signature_elem)*signatures_size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANON, 0, 0);
    if (signatures == (signature_elem*)-1) {
        fprintf(stderr, "sigref: Unable to allocate memory (%'zu bytes) for the signatures!\n", signatures_size*sizeof(signature_elem));
        exit(1);
    }
    refine_iteration++;
}

TASK_2(BDD, assign_block, BDD, sig, BDD, previous_block)
{
    assert(previous_block != mtbdd_false); // if so, incorrect call!

    // maybe do garbage collection
    sylvan_gc_test();

    if (sig == sylvan_false) {
        // slightly different handling because sylvan_false == 0
        sig = (uint64_t)-1;
    }

    // try to claim previous block number
    const uint64_t p_b = CALL(decode_block, previous_block);
    assert(p_b != 0);

    for (;;) {
        BDD cur = *(volatile BDD*)&signatures[p_b].sig;
        if (cur == sig) return previous_block;
        if (cur != 0) break;
        if (__sync_bool_compare_and_swap(&signatures[p_b].sig, 0, sig)) return previous_block;
    }

    /* claim unsuccesful, find newly added */
    uint32_t trace[SL_DEPTH];
    uint32_t loc = 0, loc_next = 0, k = SL_DEPTH-1;
    for (;;) {
        /* invariant: [loc].sig < sig */
        /* note: this is always true for loc==0 */
        signature_elem *e = &signatures[loc];
        loc_next = (*(volatile uint32_t*)&e->next[k]) & 0x7fffffff;
        if (loc_next != 0 && signatures[loc_next].sig == sig && signatures[loc_next].prev == p_b) {
            /* found */
            return CALL(encode_block, loc_next);
        } else if (loc_next != 0 && signatures[loc_next].sig == sig && signatures[loc_next].prev < p_b) {
            /* go right */
            loc = loc_next;
        } else if (loc_next != 0 && signatures[loc_next].sig < sig) {
            /* go right */
            loc = loc_next;
        } else if (k > 0) {
            /* go down */
            trace[k] = loc;
            k--;
        } else if (!(e->next[0] & 0x80000000) && __sync_bool_compare_and_swap(&e->next[0], loc_next, loc_next|0x80000000)) {
            /* locked */
            break;
        }
    }

    /* claim next item */
    const uint32_t b_nr = __sync_fetch_and_add(&next_block, 1);

    if (b_nr >= signatures_size) {
        fprintf(stderr, "Out of cheese exception, no more blocks available\n");
        exit(1);
    }

    /* fill next item */
    signature_elem *a = &signatures[b_nr];
    a->sig = sig;
    a->prev = p_b;
    a->next[0] = loc_next;
    signatures[loc].next[0] = b_nr; // TODO make the "next" pointer atomic

    /* determine height */
    uint64_t h = 1 + __builtin_clz(LACE_TRNG) / 2;
    if (h > SL_DEPTH) h = SL_DEPTH;

    /* go up and create links */
    for (k=1;k<h;k++) {
        loc = trace[k];
        for (;;) {
            signature_elem *e = &signatures[loc];
            /* note, at k>0, no locks on edges */
            uint32_t loc_next = *(volatile uint32_t*)&e->next[k];
            if (loc_next != 0 && signatures[loc_next].sig == sig && signatures[loc_next].prev < p_b) {
                loc = loc_next;
            } else if (loc_next != 0 && signatures[loc_next].sig < sig) {
                loc = loc_next;
            } else {
                a->next[k] = loc_next;
                if (__sync_bool_compare_and_swap(&e->next[k], loc_next, b_nr)) break;
            }
        }
    }

    return CALL(encode_block, b_nr);
}

TASK_4(BDD, refine_partition, BDD, dd, BDD, s_vars, BDD, ns_vars, BDD, previous_partition)
{
    // expecting dd as in s,a,B
    // expecting s_vars and ns_vars to be the set of variables of states / next states
    // expecting previous_partition as in ns,B

    if (previous_partition == sylvan_false) {
        // it had no block in the previous iteration, therefore also not now
        return sylvan_false;
    }

    if (sylvan_set_isempty(s_vars)) {
        // vars is empty, so we've reached a signature node! assign block and return
        BDD result;
        if (cache_get3(CACHE_REFINE, dd, s_vars, previous_partition, &result)) return result;
        result = CALL(assign_block, dd, previous_partition);
        cache_put3(CACHE_REFINE, dd, s_vars, previous_partition, result);
        return result;
    }

    // check if we need to do garbage collection
    sylvan_gc_test();

    /* vars != sylvan_false */
    /* dd cannot be sylvan_true - if vars != sylvan_true, then dd is in a,B */

    BDDVAR dd_var = sylvan_isconst(dd) ? 0xffffffff : sylvan_var(dd);
    BDDVAR pp_var = sylvan_var(previous_partition);
    BDDVAR s_vars_var = sylvan_set_first(s_vars);
    BDDVAR ns_vars_var = sylvan_set_first(ns_vars);

    while (s_vars_var < dd_var && ns_vars_var < pp_var) {
        s_vars = sylvan_set_next(s_vars);
        ns_vars = sylvan_set_next(ns_vars);
        if (sylvan_set_isempty(s_vars)) return CALL(refine_partition, dd, s_vars, ns_vars, previous_partition);
        s_vars_var = sylvan_set_first(s_vars);
        ns_vars_var = sylvan_set_first(ns_vars);
    }

    /* Consult cache */
    BDD result;
    if (cache_get3(CACHE_REFINE, dd, s_vars, previous_partition, &result)) {
        return result;
    }

    /* Compute cofactors */
    BDD dd_low, dd_high;
    if (s_vars_var == dd_var) {
        dd_low = sylvan_low(dd);
        dd_high = sylvan_high(dd);
    } else {
        dd_low = dd_high = dd;
    }

    BDD pp_low, pp_high;
    if (ns_vars_var == pp_var) {
        pp_low = sylvan_low(previous_partition);
        pp_high = sylvan_high(previous_partition);
    } else {
        pp_low = pp_high = previous_partition;
    }

    /* Recursive steps */
    BDD next_s_vars = sylvan_set_next(s_vars);
    BDD next_ns_vars = sylvan_set_next(ns_vars);
    bdd_refs_spawn(SPAWN(refine_partition, dd_low, next_s_vars, next_ns_vars, pp_low));
    BDD high = bdd_refs_push(CALL(refine_partition, dd_high, next_s_vars, next_ns_vars, pp_high));
    BDD low = bdd_refs_sync(SYNC(refine_partition));
    bdd_refs_pop(1);

    /* rename from s to t */
    result = sylvan_makenode(ns_vars_var, low, high);

    /* Write to cache */
    cache_put3(CACHE_REFINE, dd, s_vars, previous_partition, result);
    return result;
}

TASK_4(BDD, refine, MTBDD, signature, BDD, s_vars, BDD, ns_vars, BDD, previous_partition)
{
    prepare_refine();
    return CALL(refine_partition, signature, s_vars, ns_vars, previous_partition);
}

size_t
count_blocks()
{
    return next_block - 1;
}

void
set_signatures_size(size_t count)
{
    signatures_size = count;
}

size_t
get_next_block()
{
    return next_block++;
}

BDD
get_signature(size_t index)
{
    BDD result = signatures[index+1].sig;
    if (result == (uint64_t)-1) return sylvan_false;
    else return result;
}

void
free_refine_data()
{
    if (signatures != NULL) {
        munmap(signatures, sizeof(signature_elem)*signatures_size);
        signatures = NULL;
    }
}

TASK_IMPL_1(MTBDD, min_lts_strong, SymGame*, sym)
{
    // first get the transition relation with the priorities
    //  i.e., only keep s > uap > cap > ns
    MTBDD trans = sylvan_exists(sym->trans, sym->p_vars);
    mtbdd_refs_pushptr(&trans);

    // next, get state data, prepare blocks
    {
        block_length = sym->statebits < 25 ? sym->statebits+1 : 25; // Cap it on 2^25 : 33,554,432 blocks max
        uint32_t block_vars[block_length];
        uint32_t block_base = 1000000;    // first block BDD variables
        for (int i=0; i<block_length; i++) block_vars[i] = block_base+2*i;
        block_variables = sylvan_set_fromarray(block_vars, block_length);
        sylvan_protect(&block_variables); // perpetual
    }
    set_signatures_size(1ULL << block_length);

    // acquire ids for cache
    CACHE_ENCODE_BLOCK = cache_next_opid();
    CACHE_DECODE_BLOCK = cache_next_opid();

    // "make" variables and protect them
    MTBDD partition = mtbdd_false;
    MTBDD signature = mtbdd_false;
    mtbdd_refs_pushptr(&partition);
    mtbdd_refs_pushptr(&signature);

    // initial partition!
    partition = CALL(encode_block, get_next_block());

    // restrict to reachable states
    {
        MTBDD states = sylvan_project(sym->trans, sym->s_vars);
        mtbdd_refs_pushptr(&states);
        MTBDD s_to_ns = mtbdd_map_empty();
        mtbdd_refs_pushptr(&s_to_ns);
        MTBDD _s = sym->s_vars;
        MTBDD _n = sym->ns_vars;
        while (!mtbdd_set_isempty(_s)) {
            int sv = mtbdd_set_first(_s);
            int nv = mtbdd_set_first(_n);
            _s = mtbdd_set_next(_s);
            _n = mtbdd_set_next(_n);
            s_to_ns = mtbdd_map_add(s_to_ns, sv, sylvan_ithvar(nv));
        }
        states = sylvan_compose(states, s_to_ns);
        partition = sylvan_ite(states, partition, mtbdd_false);
        mtbdd_refs_popptr(2);
    }

    // main loop
    size_t n_blocks = 1;
    size_t iterations = 0;
    size_t old_n_blocks = 0;
    while (n_blocks != old_n_blocks) {
        old_n_blocks = n_blocks;
        iterations++;

        // std::cerr << "current partition:" << std::endl;
        // CALL(print_partition, sym, partition);

        // compute signature: s > uap > cap > B
        signature = sylvan_and_exists(trans, partition, sym->ns_vars);
        CACHE_REFINE = cache_next_opid(); // each refine needs fresh cache
        partition = CALL(refine, signature, sym->s_vars, sym->ns_vars, partition);
        n_blocks = count_blocks();
        signature = mtbdd_false; // sneaky unref
    }

    // std::cerr << "bisimulation minimisation required " << iterations << " iterations." << std::endl;
    // std::cerr << "afterwards: " << n_blocks << " blocks." << std::endl;

    mtbdd_refs_popptr(3);
    return partition;
}


uint64_t state_0_block = 0;


/**
 * Apply partition, ie assuming partition is defined on s_vars -> block_variables
 * and converts block to the ns_vars (typically short version of s_vars)
 */
TASK_4(MTBDD, apply_partition, MTBDD, dd, MTBDD, partition, MTBDD, s_vars, MTBDD, new_s_vars)
{
    // either partition on s,B and s_vars=s and b_to_s B->s
    //   or   partition on t,B and s_vars=t and b_to_s B->t
    // variables order: s < a < t < B

    // obviously, false becomes false
    if (dd == mtbdd_false) return mtbdd_false;

    // ignore states without a block
    if (partition == mtbdd_false) return mtbdd_false;

    uint32_t dd_var = dd != mtbdd_true ? sylvan_var(dd) : (uint32_t)-1;
    uint32_t partition_var = sylvan_var(partition);
    uint32_t top_var = dd_var < partition_var ? dd_var : partition_var;

    // skip in s_vars
    while (mtbdd_set_first(s_vars) < top_var) {
        s_vars = mtbdd_set_next(s_vars);
        if (mtbdd_set_isempty(s_vars)) break;
    }

    // check if we need to do garbage collection
    sylvan_gc_test();

    // no need to check s_vars / new_s_vars since we'll change CACHE_APPLY_PARTITION
    BDD result;
    if (cache_get3(CACHE_APPLY_PARTITION, dd, partition, 0, &result)) {
        return result;
    }

    // check if we just reached the end of s_vars, in which case we compute the result and return
    if (mtbdd_set_isempty(s_vars)) {
        assert(mtbdd_getvar(partition) == mtbdd_set_first(block_variables));

        uint64_t block_no = CALL(decode_block, partition);
        // std::cerr << "pre: block " << block_no << std::endl;
        if (block_no == state_0_block) block_no = 1;
        else if (block_no == 1) block_no = state_0_block;
        block_no--; // block start with block 1...

        // HERE modify for state "0"... ensure state "0" is block "0"!

        // std::cerr << "got block " << block_no << " and dd " << dd << std::endl;

        result = mtbdd_true;
        mtbdd_refs_pushptr(&result);

        int count = mtbdd_set_count(new_s_vars);
        for (int j=0; j<count; j++) {
            // assume consecutive variables...
            int var = mtbdd_set_first(new_s_vars)+count-j-1;

            if (block_no & (1ULL<<j)) {
                result = mtbdd_makenode(var, mtbdd_false, result);
            } else {
                result = mtbdd_makenode(var, result, mtbdd_false);
            }
        }

        result = sylvan_and(result, dd);
        mtbdd_refs_popptr(1);

        /* cache and return */
        cache_put3(CACHE_APPLY_PARTITION, dd, partition, 0, result);
        return result;
    }

    // compute cofactors
    BDD dd_low, dd_high;
    if (dd_var == top_var) {
        dd_low = sylvan_low(dd);
        dd_high = sylvan_high(dd);
    } else {
        dd_low = dd_high = dd;
    }

    BDD partition_low, partition_high;
    if (partition_var == top_var) {
        partition_low = sylvan_low(partition);
        partition_high = sylvan_high(partition);
    } else {
        partition_low = partition_high = partition;
    }

    // compute recursive results
    mtbdd_refs_spawn(SPAWN(apply_partition, dd_low, partition_low, s_vars, new_s_vars));
    BDD high = mtbdd_refs_push(CALL(apply_partition, dd_high, partition_high, s_vars, new_s_vars));
    BDD low = mtbdd_refs_push(mtbdd_refs_sync(SYNC(apply_partition)));

    // Now depending on whether the current var is in s_vars or not we return...
    bool is_s_var = mtbdd_set_first(s_vars) == top_var;
    if (is_s_var) {
        // in this case either
        //  - s > a > B
        //    then the recursive results are full s>a>B that need to be ORed
        //  - .. > a > s > B
        //    then again the recursive results are full states that need to be ORed
        result = sylvan_or(low, high);
    } else {
        // in this case either
        //  - s > a > B
        //    then we should not even be here...
        //  - .. > a > s > B
        //    then we are /before/ s and can just create the node
        assert(mtbdd_isleaf(low) || top_var < mtbdd_getvar(low));
        assert(mtbdd_isleaf(high) || top_var < mtbdd_getvar(high));
        result = sylvan_makenode(top_var, low, high);
    }

    mtbdd_refs_pop(2);

    /* cache result */
    cache_put3(CACHE_APPLY_PARTITION, dd, partition, 0, result);

    return result;
}

VOID_TASK_IMPL_3(minimize, SymGame*, sym, MTBDD, partition, bool, verbose)
{
    // first, find out where state 0 goes
    {
        MTBDD _p = partition;
        while (mtbdd_getvar(_p) < mtbdd_set_first(block_variables)) _p = mtbdd_getlow(_p);
        state_0_block = CALL(decode_block, _p);
    }

    // get number of blocks
    size_t no_blocks = count_blocks();

    // compute how many ns variables are needed
    int no_ns_vars = 1;
    while (no_blocks > (1ULL << no_ns_vars)) no_ns_vars++;

    if (verbose) {
        std::cerr << "after bisimulation minimisation: " << no_blocks << " blocks." << std::endl;
        // std::cerr << "state 0 is in block " << state_0_block << "; there are " << no_blocks << " blocks." << std::endl;
        // std::cerr << "now using " << no_ns_vars << " state variable(s) instead of " << sym->statebits << "." << std::endl;
    }

    MTBDD block_s_vars = mtbdd_set_empty();
    MTBDD block_ns_vars = mtbdd_set_empty();
    MTBDD ns_to_s_map = mtbdd_map_empty();

    mtbdd_refs_pushptr(&block_s_vars);
    mtbdd_refs_pushptr(&block_ns_vars);
    mtbdd_refs_pushptr(&ns_to_s_map);
    mtbdd_refs_pushptr(&partition);

    {
        // Construct new s_vars
        uint32_t first = mtbdd_set_first(sym->s_vars);
        for (int i=0; i<no_ns_vars; i++) block_s_vars = mtbdd_set_add(block_s_vars, first+i);
    }

    {
        // Construct new ns_vars
        uint32_t first = mtbdd_set_first(sym->ns_vars);
        for (int i=0; i<no_ns_vars; i++) block_ns_vars = mtbdd_set_add(block_ns_vars, first+i);
    }

    {
        // Construct ns_to_s (on original state/next_state variables)
        MTBDD _s = sym->s_vars;
        MTBDD _n = sym->ns_vars;
        while (!mtbdd_set_isempty(_s)) {
            int sv = mtbdd_set_first(_s);
            int nv = mtbdd_set_first(_n);
            _s = mtbdd_set_next(_s);
            _n = mtbdd_set_next(_n);
            ns_to_s_map = mtbdd_map_add(ns_to_s_map, nv, sylvan_ithvar(sv));
        }
    }

    // run on next state (only trans)
    CACHE_APPLY_PARTITION = cache_next_opid();
    // if (verbose) std::cerr << "applying partition to next states of transition relation" << std::endl;
    sym->trans = CALL(apply_partition, sym->trans, partition, sym->ns_vars, block_ns_vars);
    // sym->print_trans();
    
    // change partition to s -> B
    // if (verbose) std::cerr << "preparing partition on current states" << std::endl;
    partition = mtbdd_compose(partition, ns_to_s_map);

    // run on state (trans and strategies)
    CACHE_APPLY_PARTITION = cache_next_opid();
    // if (verbose) std::cerr << "applying partition to current states of transition relation" << std::endl;
    sym->trans = CALL(apply_partition, sym->trans, partition, sym->s_vars, block_s_vars);
    // sym->print_trans();
    // if (verbose) std::cerr << "applying partition to states of strategies" << std::endl;
    sym->strategies = CALL(apply_partition, sym->strategies, partition, sym->s_vars, block_s_vars);
    // sym->print_strategies();

    sym->s_vars = block_s_vars;
    sym->ns_vars = block_ns_vars;
    sym->pns_vars = mtbdd_set_addall(sym->p_vars, sym->ns_vars);
    sym->statebits = no_ns_vars;
    
    mtbdd_refs_popptr(4);
}


VOID_TASK_IMPL_2(print_partition, SymGame*, game, MTBDD, partition)
{
    uint8_t arr[game->statebits+1];
    MTBDD lf = mtbdd_enum_all_first(partition, game->ns_vars, arr, NULL);
    while (lf != mtbdd_false) {
        int idx=0;
        // decode state
        int state = 0;
        for (int i=0; i<game->statebits; i++) {
            state <<= 1;
            if (arr[idx++]) state |= 1;
        }
        // std::cerr << "top var of lf is " << mtbdd_getvar(lf) << std::endl;
        assert(mtbdd_getvar(lf) == mtbdd_set_first(block_variables));
        // decode block
        uint64_t block = CALL(decode_block, lf);
        std::cerr << "state " << state << ": " << block << std::endl;
        lf = mtbdd_enum_all_next(partition, game->ns_vars, arr, NULL);
    }
}



