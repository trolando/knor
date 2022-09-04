/**************************************************************************
 * Copyright Tom van Dijk
 *************************************************************************/

#include <aigmaker.hpp>
#include <knor.hpp>
#include <set>
#include <cstring> // for strcpy

using namespace sylvan;


static void
collect_inter(MTBDD bdd, uint32_t firstvar, std::set<MTBDD> &res)
{
    // TODO needs caching (visited or not)
    if (bdd == mtbdd_false) {
        return;
    }
    else if (mtbdd_isleaf(bdd)) {
        res.insert(bdd);
    } else {
        if (mtbdd_getvar(bdd) < firstvar) {
            collect_inter(mtbdd_gethigh(bdd), firstvar, res);
            collect_inter(mtbdd_getlow(bdd), firstvar, res);
        } else {
            res.insert(bdd);
        }
    }
}


TASK_3(MTBDD, collect_ending, MTBDD, bdd, uint32_t, firstvar, MTBDD, leaf)
{
    // TODO needs caching
    if (bdd == leaf) {
        return mtbdd_true;
    } else if (mtbdd_isleaf(bdd)) {
        return mtbdd_false;
    } else if (mtbdd_getvar(bdd) >= firstvar) {
        return mtbdd_false;
    } else {
        uint32_t var = mtbdd_getvar(bdd);
        mtbdd_refs_spawn(SPAWN(collect_ending, mtbdd_gethigh(bdd), firstvar, leaf));
        MTBDD low = mtbdd_refs_push(CALL(collect_ending, mtbdd_getlow(bdd), firstvar, leaf));
        MTBDD high = mtbdd_refs_push(mtbdd_refs_sync(SYNC(collect_ending)));
        MTBDD res = mtbdd_makenode(var, low, high);
        mtbdd_refs_pop(2);
        return res;
    }
}


AIGmaker::AIGmaker(HoaData *data, SymGame *game) : data(data), game(game)
{
    a = NULL;
}

AIGmaker::~AIGmaker()
{
    if (a != NULL) {
        delete[] uap_to_lit;
        delete[] caps;
    }
}

int
AIGmaker::makeand(int rhs0, int rhs1)
{
    if (rhs1 < rhs0) {
        int tmp = rhs0;
        rhs0 = rhs1;
        rhs1 = tmp;
    }

    if (rhs0 == 0) return 0;
    if (rhs0 == 1) return rhs1;

    uint64_t cache_key = rhs1;
    cache_key <<= 32;
    cache_key |= rhs0;
    auto c = cache.find(cache_key);
    if (1 and c != cache.end()) {
        return c->second;
    } else {
        aiger_add_and(a, lit, rhs0, rhs1);
        cache[cache_key] = lit;
        lit += 2;
        return lit-2;
    }
}

void
AIGmaker::simplify_and(std::deque<int> &gates)
{
    // for each pair of gates in gates, check the cache
    for (auto first = gates.begin(); first != gates.end(); ++first) {
        for (auto second = first + 1; second != gates.end(); ++second) {
            int left = *first;
            int right = *second;
            if (left > right) std::swap(left, right);
            uint64_t cache_key = right;
            cache_key <<= 32;
            cache_key |= left;
            auto c = cache.find(cache_key);
            if (c != cache.end()) {
                //gates.erase(std::remove_if(gates.begin(), gates.end(), [=](int x){return x==left or x==right;}),
                //        gates.end());
                gates.erase(second);
                gates.erase(first);
                gates.push_back(c->second);
                simplify_and(gates);
                return;
            }
        }
    }
}

void
AIGmaker::simplify_or(std::deque<int> &gates)
{
    // for each pair of gates in gates, check the cache
    for (auto first = gates.begin(); first != gates.end(); ++first) {
        for (auto second = first + 1; second != gates.end(); ++second) {
            int left = aiger_not(*first);
            int right = aiger_not(*second);
            if (left > right) std::swap(left, right);
            uint64_t cache_key = right;
            cache_key <<= 32;
            cache_key |= left;
            auto c = cache.find(cache_key);
            if (c != cache.end()) {
                gates.erase(second);
                gates.erase(first);
                gates.push_back(aiger_not(c->second));
                simplify_or(gates);
                return;
            }
        }
    }
}

int
AIGmaker::bdd_to_aig_isop(MTBDD bdd)
{
    if (verbose) {
        // std::cerr << "running isop for BDD with " << mtbdd_nodecount(bdd) << " nodes." << std::endl;
    }
    MTBDD bddres;
    ZDD isop = zdd_isop(bdd, bdd, &bddres);
    zdd_protect(&isop);
    // no need to reference the result...
    assert(bdd == bddres);
    assert(bdd == zdd_cover_to_bdd(isop));
    if (verbose) {
        // std::cerr << "isop has " << (long)zdd_pathcount(isop) << " terms and " << zdd_nodecount(&isop, 1) << " nodes." << std::endl;
    }

    int res = bdd_to_aig_cover(isop);
    zdd_unprotect(&isop);
    return res;
}

int
AIGmaker::make_and(std::deque<int> &gates)
{
    while (!gates.empty()) {
        int last = gates.front();
        gates.pop_front();
        if (!gates.empty()) {
            int last2 = gates.front();
            gates.pop_front();
            int new_gate = makeand(last, last2);
            gates.push_back(new_gate);
        } else {
            return last;
        }
    }
    return aiger_false;
}

int
AIGmaker::make_or(std::deque<int> &gates)
{
    while (!gates.empty()) {
        int last = gates.front();
        gates.pop_front();
        if (!gates.empty()) {
            int last2 = gates.front();
            gates.pop_front();
            int new_gate = aiger_not(makeand(aiger_not(last), aiger_not(last2)));
            gates.push_back(new_gate);
        } else {
            return last;
        }
    }
    return aiger_false;
}

int
AIGmaker::bdd_to_aig_cover_sop(ZDD cover)
{
    if (cover == zdd_true) return aiger_true;
    if (cover == zdd_false) return aiger_false;

    // a product could consist of all variables, and a -1 to denote
    //  the end of the product
    int product[game->statebits+game->uap_count+1] = { 0 };

    // a queue that stores all products, which will need to be summed
    std::deque<int> products;

    ZDD res = zdd_cover_enum_first(cover, product);
    while (res != zdd_false) {
        //  containing subproducts in the form of gates
        std::deque<int> gates;

        for (int i=0; product[i] != -1; i++) {
            int the_lit = var_to_lit[product[i]/2];
            if (product[i]&1) the_lit = aiger_not(the_lit);
            gates.push_back(the_lit);
        }

        // simplify_and(gates);

        // while we still have subproducts we need to AND together
        while (!gates.empty()) {
            int last = gates.front();
            gates.pop_front();
            if (!gates.empty()) {
                int last2 = gates.front();
                gates.pop_front();
                int new_gate = makeand(last, last2);
                gates.push_back(new_gate);
            } else {
                products.push_back(last);
            }
        }
        res = zdd_cover_enum_next(cover, product); // go to the next product
    }

    // products queue should now be full of complete products that need to be summed

    // simplify_or(products);

    while (!products.empty()) {
        int product1 = products.front();
        products.pop_front();
        if (!products.empty()) {
            int product2 = products.front();
            products.pop_front();
            int summed_product = aiger_not(makeand(aiger_not(product1), aiger_not(product2)));
            products.push_back(summed_product);
        } else { // product1 is the final sum of all products
            return product1;
        }
    }

    return aiger_false; // should be unreachable, tbh.
}

int
AIGmaker::bdd_to_aig_cover(ZDD cover)
{
    if (cover == zdd_true) return aiger_true;
    if (cover == zdd_false) return aiger_false;

    auto it = mapping.find(cover);
    if (it != mapping.end()) {
        return it->second;
    }

    int the_var = zdd_getvar(cover);
    int the_lit = var_to_lit[the_var/2];
    if (the_var & 1) the_lit = aiger_not(the_lit);

    ZDD low = zdd_getlow(cover);
    ZDD high = zdd_gethigh(cover);

    int res = the_lit;

    if (high != zdd_true) {
        auto x = bdd_to_aig_cover(high);
        res = makeand(res, x);
    }

    if (low != zdd_false) {
        auto x = bdd_to_aig_cover(low);
        res = aiger_not(makeand(aiger_not(res), aiger_not(x)));
    }

    mapping[cover] = res;
    return res;
}

int
AIGmaker::bdd_to_aig(MTBDD bdd)
{
    if (bdd == mtbdd_true) return aiger_true;
    if (bdd == mtbdd_false) return aiger_false;

    bool comp = false;
    if (bdd & sylvan_complement) {
        bdd ^= sylvan_complement;
        comp = true;
    }

    auto it = mapping.find(bdd);
    if (it != mapping.end()) {
        return comp ? aiger_not(it->second) : it->second;
    }

    int the_lit = var_to_lit[mtbdd_getvar(bdd)];

    MTBDD low = mtbdd_getlow(bdd);
    MTBDD high = mtbdd_gethigh(bdd);

    int res;

    if (low == mtbdd_false) {
        // only high (value 1)
        if (high == mtbdd_true) {
            // actually this is the end, just the lit
            res = the_lit;
        } else {
            // AND(the_lit, ...)
            int rhs0 = the_lit;
            int rhs1 = bdd_to_aig(high);
            res = makeand(rhs0, rhs1);
        }
    } else if (high == mtbdd_false) {
        // only low (value 0)
        if (low == mtbdd_true) {
            // actually this is the end, just the lit, negated
            res = aiger_not(the_lit);
        } else {
            // AND(not the_lit, ...)
            int rhs0 = aiger_not(the_lit);
            int rhs1 = bdd_to_aig(low);
            res = makeand(rhs0, rhs1);
        }
    } else {
        // OR(low, high) == ~AND(~AND(the_lit, ...), ~AND(~the_lit, ...))
        int lowres = bdd_to_aig(low);
        int highres = bdd_to_aig(high);
        int rhs0 = aiger_not(makeand(aiger_not(the_lit), lowres));
        int rhs1 = aiger_not(makeand(the_lit, highres));
        res = aiger_not(makeand(rhs0, rhs1));
    }

    mapping[bdd] = res;

    return comp ? aiger_not(res) : res;
}

void
AIGmaker::process()
{
    a = aiger_init();
    lit = 2;

    // Set which APs are controllable in the bitset controllable
    pg::bitset controllable(data->noAPs);
    for (int i=0; i<data->noCntAPs; i++) {
        controllable[data->cntAPs[i]] = 1;
    }

    uap_to_lit = new int[game->uap_count];
    caps = new char*[game->cap_count];

    int uap_idx = 0;
    int cap_idx = 0;
    for (int i=0; i<data->noAPs; i++) {
        if (!controllable[i]) {
            uap_to_lit[uap_idx] = lit;
            aiger_add_input(a, lit, data->aps[i]);
            uap_idx++;
            lit += 2;
        } else {
            caps[cap_idx++] = data->aps[i];
        }
    }

    /* make var_to_lit for uncontrollable APs */
    {
        MTBDD _vars = game->uap_vars;
        for (int i=0; i<game->uap_count; i++) {
            uint32_t bddvar = mtbdd_set_first(_vars);
            _vars = mtbdd_set_next(_vars);
            var_to_lit[bddvar] = uap_to_lit[i];
        }
    }

    if (onehot) {
        MTBDD states = sylvan_project(game->trans, game->s_vars);
        mtbdd_protect(&states);

        /* prepare state_to_lit */
        {
            uint8_t state_arr[game->statebits];
            MTBDD lf = mtbdd_enum_all_first(states, game->s_vars, state_arr, NULL);
            while (lf != mtbdd_false) {
                // decode state
                int state = 0;
                for (int i=0; i<game->statebits; i++) {
                    state <<= 1;
                    if (state_arr[i]) state |= 1;
                }
                // give the state a literal
                // std::cerr << "state " << state << " gets latch literal " << lit << std::endl;
                state_to_lit[state] = lit;
                lit += 2;
                // next state
                lf = mtbdd_enum_all_next(states, game->s_vars, state_arr, NULL);
            }
        }

        // do each controllable AP (output signals)
        MTBDD cap_bdd = mtbdd_false;
        mtbdd_protect(&cap_bdd);

        MTBDD s = mtbdd_false;
        mtbdd_protect(&s);

        for (int i=0; i<game->cap_count; i++) {
            cap_bdd = sylvan_ithvar(mtbdd_set_first(game->cap_vars)+i);
            // keep just s and u... get rid of other cap variables
            cap_bdd = sylvan_and_exists(game->strategies, cap_bdd, game->cap_vars);

            // std::cerr << "for controllable ap " << i << std::endl;

            // this gives source states and UAP for this cap
            std::set<MTBDD> uaps;
            collect_inter(cap_bdd, mtbdd_set_first(game->uap_vars), uaps);

            std::deque<int> terms;

            for (MTBDD uap : uaps) {
                s = RUN(collect_ending, cap_bdd, mtbdd_set_first(game->uap_vars), uap);

                std::vector<int> source_states;
                std::deque<int> source_gates;

                uint8_t state_arr_2[game->statebits];
                MTBDD lf2 = mtbdd_enum_all_first(s, game->s_vars, state_arr_2, NULL);
                while (lf2 != mtbdd_false) {
                    // decode state
                    int from = 0;
                    for (int i=0; i<game->statebits; i++) {
                        from <<= 1;
                        if (state_arr_2[i]) from |= 1;
                    }
                    source_states.push_back(from);
                    // if state 0, take inverse
                    source_gates.push_back(from == 0 ? aiger_not(state_to_lit[from]) : state_to_lit[from]);
                    lf2 = mtbdd_enum_all_next(s, game->s_vars, state_arr_2, NULL);
                }

                //for (int from : source_states) {
                    //std::cerr << "source state " << from << " with UAP " << ((uap&sylvan_complement)?"~":"") << (uap&~sylvan_complement) << std::endl;
                //}

                int aig_uap = isop ? bdd_to_aig_isop(uap) : bdd_to_aig(uap);
                int aig_states = make_or(source_gates);
                terms.push_back(makeand(aig_uap, aig_states));
            }

            int result = make_or(terms);
            aiger_add_output(a, result, caps[i]); // simple, really
        }

        MTBDD su_vars = mtbdd_set_addall(game->p_vars, game->cap_vars);
        mtbdd_protect(&su_vars);

        // full is: s > u > ns
        MTBDD full = mtbdd_and_exists(game->strategies, game->trans, su_vars);
        mtbdd_protect(&full);

        // long noStates = (long)sylvan_satcount(states, game->s_vars);
        // std::cerr << "number of states: " << noStates << std::endl;

        uint8_t state_arr[game->statebits];
        MTBDD lf = mtbdd_enum_all_first(states, game->s_vars, state_arr, NULL);
        while (lf != mtbdd_false) {
            // decode state
            int state = 0;
            for (int i=0; i<game->statebits; i++) {
                state <<= 1;
                if (state_arr[i]) state |= 1;
            }

            // encode state using NS vars
            cap_bdd = encode_state(state, game->statebits, mtbdd_set_first(game->ns_vars));
            // keep s > u of this state
            cap_bdd = sylvan_and_exists(full, cap_bdd, game->ns_vars);

            // std::cerr << "state " << state << /*" to has " << s <<*/ std::endl;

            std::set<MTBDD> uaps;
            collect_inter(cap_bdd, mtbdd_set_first(game->uap_vars), uaps);

            std::deque<int> terms;
            for (MTBDD uap : uaps) {
                s = RUN(collect_ending, cap_bdd, mtbdd_set_first(game->uap_vars), uap);

                std::vector<int> source_states;
                std::deque<int> source_gates;

                uint8_t state_arr_2[game->statebits];
                MTBDD lf2 = mtbdd_enum_all_first(s, game->s_vars, state_arr_2, NULL);
                while (lf2 != mtbdd_false) {
                    // decode state
                    int from = 0;
                    for (int i=0; i<game->statebits; i++) {
                        from <<= 1;
                        if (state_arr_2[i]) from |= 1;
                    }
                    source_states.push_back(from);
                    // if state 0, take inverse
                    source_gates.push_back(from == 0 ? aiger_not(state_to_lit[from]) : state_to_lit[from]);
                    lf2 = mtbdd_enum_all_next(s, game->s_vars, state_arr_2, NULL);
                }

                //for (int from : source_states) {
                    // std::cerr << "source state " << from << " with UAP " << ((uap&sylvan_complement)?"~":"") << (uap&~sylvan_complement) << std::endl;
                //}

                int aig_uap = isop ? bdd_to_aig_isop(uap) : bdd_to_aig(uap);
                int aig_states = make_or(source_gates);
                int result = makeand(aig_uap, aig_states);
                terms.push_back(result);
            }

            int result = make_or(terms);
            if (state == 0) result = aiger_not(result);
            aiger_add_latch(a, state_to_lit[state], result, "");
            lf = mtbdd_enum_all_next(states, game->s_vars, state_arr, NULL);
        }

        mtbdd_unprotect(&states);
        mtbdd_unprotect(&cap_bdd);
        mtbdd_unprotect(&s);
        mtbdd_unprotect(&full);
        mtbdd_unprotect(&su_vars);
    } else {
        // NOT onehot

        /* make state_to_lit for latches */
        for (int i=0; i<game->statebits; i++) {
            state_to_lit[i] = lit;
            lit += 2;
        }

        /* make var_to_lit for states */
        MTBDD _vars = game->s_vars;
        for (int i=0; i<game->statebits; i++) {
            uint32_t bddvar = mtbdd_set_first(_vars);
            _vars = mtbdd_set_next(_vars);
            var_to_lit[bddvar] = state_to_lit[i];
        }

        MTBDD* cap_bdds;   // contains the solution: controllable ap bdds: state -> uap -> B
        MTBDD* state_bdds; // contains the solution: state bit bdds      : state -> uap -> B

        // compute bdds for the controllable APs
        cap_bdds = new MTBDD[game->cap_count];
        for (int i=0; i<game->cap_count; i++) {
            cap_bdds[i] = sylvan_ithvar(mtbdd_set_first(game->cap_vars)+i);
            mtbdd_protect(&cap_bdds[i]);
            // keep just s and u... get rid of other cap variables
            cap_bdds[i] = sylvan_and_exists(game->strategies, cap_bdds[i], game->cap_vars);
        }

        // su_vars is priority and cap, which is to be removed...
        MTBDD su_vars = mtbdd_set_addall(game->p_vars, game->cap_vars);
        mtbdd_protect(&su_vars);

        // full is: s > u > ns
        MTBDD full = mtbdd_and_exists(game->strategies, game->trans, su_vars);
        mtbdd_protect(&full);

        // compute bdds for each state variable
        state_bdds = new MTBDD[game->statebits];
        for (int i=0; i<game->statebits; i++) {
            state_bdds[i] = sylvan_ithvar(mtbdd_set_first(game->ns_vars) + i);
            mtbdd_protect(&state_bdds[i]);
            // keep just s and u... get rid of other cap variables
            state_bdds[i] = sylvan_and_exists(full, state_bdds[i], game->ns_vars);
        }

        mtbdd_unprotect(&full);
        mtbdd_unprotect(&su_vars);

        // if ISOP, first convert all cap bdds etc to covers
        if (isop) {
            ZDD* cap_zdds = new ZDD[game->cap_count];
            for (int i=0; i<game->cap_count; i++) {
                cap_zdds[i] = zdd_false;
                zdd_protect(&cap_zdds[i]);

                MTBDD bddres;
                cap_zdds[i] = zdd_isop(cap_bdds[i], cap_bdds[i], &bddres);
                assert(bddres == cap_bdds[i]);
                if (verbose) {
                    std::cerr << "isop has " << (long)zdd_pathcount(cap_zdds[i]) << " terms and " << zdd_nodecount(&cap_zdds[i], 1) << " nodes." << std::endl;
                }
            }

            ZDD* state_zdds = new ZDD[game->statebits];
            for (int i=0; i<game->statebits; i++) {
                state_zdds[i] = zdd_false;
                zdd_protect(&state_zdds[i]);

                MTBDD bddres;
                state_zdds[i] = zdd_isop(state_bdds[i], state_bdds[i], &bddres);
                assert(bddres == state_bdds[i]);
                if (verbose) {
                    std::cerr << "isop has " << (long)zdd_pathcount(state_zdds[i]) << " terms and " << zdd_nodecount(&state_zdds[i], 1) << " nodes." << std::endl;
                }
            }

            for (int i=0; i<game->cap_count; i++) {
                int res = bdd_to_aig_cover(cap_zdds[i]);
                aiger_add_output(a, res, caps[i]); // simple, really
            }
            for (int i=0; i<game->statebits; i++) {
                int res = bdd_to_aig_cover(state_zdds[i]);
                aiger_add_latch(a, state_to_lit[i], res, "");
            }
        } else {
            for (int i=0; i<game->cap_count; i++) {
                int res = isop ? bdd_to_aig_isop(cap_bdds[i]) : bdd_to_aig(cap_bdds[i]);
                aiger_add_output(a, res, caps[i]); // simple, really
            }
            for (int i=0; i<game->statebits; i++) {
                int res = isop ? bdd_to_aig_isop(state_bdds[i]) : bdd_to_aig(state_bdds[i]);
                aiger_add_latch(a, state_to_lit[i], res, "");
            }
        }

        for (int i=0; i<game->cap_count; i++) mtbdd_unprotect(&cap_bdds[i]);
        for (int i=0; i<game->statebits; i++) mtbdd_unprotect(&state_bdds[i]);
        delete[] cap_bdds;
        delete[] state_bdds;
    }
}

void
AIGmaker::writeAscii(FILE* out)
{
    aiger_write_to_file(a, aiger_ascii_mode, out);
}

void
AIGmaker::writeBinary(FILE* out)
{
    aiger_write_to_file(a, aiger_binary_mode, out);
}

// commands taken from 'alias compress2rs' from 'abc.rc' file
const std::vector<std::string> AIGmaker::compressCommands ({
    "balance -l",
    "resub -K 6 -l",
    "rewrite -l",
    "resub -K 6 -N 2",
    "refactor -l",
    "resub -K 8 -l",
    "balance -l",
    "resub -K 8 -N 2 -l",
    "rewrite -l",
    "resub -K 10 -l",
    "rewrite -z -l",
    "resub -K 10 -N 2 -l",
    "balance -l",
    "resub -K 12 -l",
    "refactor -z -l",
    "resub -K 12 -N 2 -l",
    "balance -l",
    "rewrite -z -l",
    "balance -l",
    "dc2"
});

// commands taken from 'alias compress2' from 'abc.rc' file
/*
const std::vector<std::string> AIGmaker::compressCommands ({
    "balance -l",
    "rewrite -l",
    "refactor -l",
    "balance -l",
    "rewrite -l",
    "rewrite -z -l",
    "balance -l",
    "refactor -z -l",
    "rewrite -z -l",
    "balance -l"
});
*/

void
AIGmaker::compress()
{
    Abc_Start();
    Abc_Frame_t* pAbc = Abc_FrameGetGlobalFrame();

    writeToAbc(pAbc);

    // compress until convergence
    int new_num_nodes = getAbcNetworkSize(pAbc);
    int old_num_nodes = new_num_nodes + 1;
    while (new_num_nodes > 0 && new_num_nodes < old_num_nodes) {
        executeCompressCommands(pAbc);
        old_num_nodes = new_num_nodes;
        new_num_nodes = getAbcNetworkSize(pAbc);
        // std::cerr << "nodes after compress run: " << new_num_nodes << std::endl;
        if ((old_num_nodes-new_num_nodes)<old_num_nodes/40) break; // 2.5% improvement or better pls
    }

    readFromAbc(pAbc);

    Abc_Stop();
}

void
AIGmaker::drewrite()
{
    Abc_Start();
    Abc_Frame_t* pAbc = Abc_FrameGetGlobalFrame();

    writeToAbc(pAbc);

    // compress until convergence
    int new_num_nodes = getAbcNetworkSize(pAbc);
    int old_num_nodes = new_num_nodes + 1;
    while (new_num_nodes > 0 && new_num_nodes < old_num_nodes) {
        executeAbcCommand(pAbc, "drw");
        executeAbcCommand(pAbc, "balance");
        executeAbcCommand(pAbc, "drf");
        executeAbcCommand(pAbc, "dc2");
        old_num_nodes = new_num_nodes;
        new_num_nodes = getAbcNetworkSize(pAbc);
        // std::cerr << "nodes after compress run: " << new_num_nodes << std::endl;
        if ((old_num_nodes-new_num_nodes)<old_num_nodes/100) break; // 1% improvement or better pls
    }

    readFromAbc(pAbc);

    Abc_Stop();
}

void AIGmaker::executeAbcCommand(Abc_Frame_t* pAbc, const std::string command) const {
    if (Cmd_CommandExecute( pAbc, command.c_str())) {
        throw std::runtime_error("Cannot execute ABC command: " + command);
    }
    if (verbose) std::cerr << "after " << command << ": " << getAbcNetworkSize(pAbc) << std::endl;
}

void AIGmaker::executeCompressCommands(Abc_Frame_t* pAbc) const {
    for (const auto& command : compressCommands) {
        executeAbcCommand(pAbc, command);
    }
}

int AIGmaker::getAbcNetworkSize(Abc_Frame_t* pAbc) const {
    Abc_Ntk_t* pNtk = Abc_FrameReadNtk(pAbc);
    return Abc_NtkNodeNum(pNtk);
}

int AIGmaker::getTmpFile(char* tmp_filename) const {
    std::strcpy(tmp_filename, "knor.XXXXXX");
    int fd = mkstemp(tmp_filename);
    if (fd == -1) {
        throw std::runtime_error("Could not create temporary file: " + std::string(tmp_filename));
    }
    return fd;
}

void AIGmaker::writeToAbc(Abc_Frame_t* pAbc) const {
    char tmp_filename[256];
    int fd = getTmpFile(tmp_filename);

    // write AIGER out to be read by ABC
    FILE* file = fdopen(fd, "w");
    if (file == nullptr) {
        throw std::runtime_error("Could not open temporary file: " + std::string(tmp_filename));
    }
    int write_result = aiger_write_to_file(a, aiger_binary_mode, file);
    fclose(file);
    if (write_result == 0) {
        throw std::runtime_error("Could not write AIGER circuit to file: " + std::string(tmp_filename));
    }

    std::stringstream read_command;
    read_command << "read_aiger " << tmp_filename;
    executeAbcCommand(pAbc, read_command.str());

    std::remove(tmp_filename);
}

void AIGmaker::readFromAbc(Abc_Frame_t* pAbc) {
    char tmp_filename[256];
    int fd = getTmpFile(tmp_filename);

    std::stringstream write_command;
    write_command << "write_aiger -s " << tmp_filename;
    executeAbcCommand(pAbc, write_command.str());

    // read AIGER back, delete comments added by ABC
    FILE* file = fdopen(fd, "r");
    if (file == nullptr) {
        throw std::runtime_error("Could not open temporary file: " + std::string(tmp_filename));
    }
    // read_aiger
    aiger_reset(a);
    a = aiger_init();
    const char* read_result = aiger_read_from_file(a, file);
    fclose(file);
    std::remove(tmp_filename);
    if (read_result != nullptr) {
        throw std::runtime_error("Could not read AIGER circuit from file: " + std::string(tmp_filename) + ": " + std::string(read_result));
    }
    aiger_delete_comments(a);
}

