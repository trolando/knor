/**************************************************************************
 * Copyright Tom van Dijk
 *************************************************************************/

#include <abcminimization.hpp>
#include <aigmaker.hpp>
#include <bddtools.hpp>
#include <knor.hpp>
#include <set>
#include <cstring> // for strcpy

// FIXME the caches should be wiped on gc, they are redundant anyway
// the important bit happens in makeand to ensure no double AIGs are made...
#define CACHES 0  // if set to 0, this will NOT use a cache on makeand!
// note that the "sop" method does not require caches anyway

using namespace sylvan;

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
    if (c != cache.end()) {
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

    int res = bdd_to_aig_cover_sop(isop);
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
            int the_lit = var_to_lit.at(product[i]/2);
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
    if (CACHES && it != mapping.end()) {
        return it->second;
    }

    int the_var = zdd_getvar(cover);
    int the_lit = var_to_lit.at(the_var/2);
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
    if (CACHES && it != mapping.end()) {
        return comp ? aiger_not(it->second) : it->second;
    }

    auto the_lit = var_to_lit.at(mtbdd_getvar(bdd));

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
AIGmaker::reduce(std::vector<std::vector<int>>& system, bool is_or)
{
    std::map<int, int> pop;
    for (auto& x : system) {
        if (x.size() == 1) continue;
        for (auto& y : x) {
            pop[y] += 1;
        }
    }

    std::map<int, int> pop2;

    while (true) {
        /*
        for (auto& x : system) {
            for (auto& y : x) std::cerr << y << " ";
            std::cerr << std::endl << std::flush;
        }
        for (auto& x : pop) {
            std::cerr << x.first << ": " << x.second << std::endl;
        }
        */

        if (pop.empty()) break;
        auto best = std::max_element(pop.begin(), pop.end(),
                 [] (const auto& a, const auto& b) -> bool { return a.second < b.second; } );
        auto first = best->first;

        pop2.clear();
        for (auto& x : system) {
            if (std::find(x.begin(), x.end(), first) != x.end()) {
                for (auto& y : x) {
                    if (y != first) pop2[y] += 1;
                }
            }
        }
        best = std::max_element(pop2.begin(), pop2.end(),
                 [] (const auto& a, const auto& b) -> bool { return a.second < b.second; } );
        auto second = best->first;

        // create new gate
        int m = 0, n = 0;
        auto gate = is_or ? aiger_not(makeand(aiger_not(first), aiger_not(second))) : makeand(first, second);
        for (auto& x : system) {
            auto it_f = std::find(x.begin(), x.end(), first);
            if (it_f != x.end()) {
                auto it_s = std::find(x.begin(), x.end(), second);
                if (it_s != x.end()) {
                    x.erase(std::remove_if(x.begin(), x.end(), [&] (const auto& a) -> bool { return a == first or a == second; }), x.end());
                    x.erase(std::remove(x.begin(), x.end(), second), x.end());
                    x.push_back(gate);
                    m++;
                    if (x.size() > 1) n++;
                }
            }
        }
        if ((pop[first] -= m) == 0) pop.erase(first);
        if ((pop[second] -= m) == 0) pop.erase(second);
        if (n > 0) pop[gate] += n;
    }
}

void
AIGmaker::process_sop()
{
    // initialize aiger
    a = aiger_init();

    // we start with literal 2 (because 0/1 is reserved)
    lit = 2;

    // set which APs are controllable in the bitset controllable
    pg::bitset controllable(data->noAPs);
    for (int i=0; i<data->noCntAPs; i++) {
        controllable[data->cntAPs[i]] = 1;
    }

    // initialize array of input signals to literal
    uap_to_lit = new int[game->uap_count];

    // initialize array of output signal labels
    caps = new char*[game->cap_count];

    // fill the two arrays
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

    // make var_to_lit for uncontrollable APs
    {
        MTBDD _vars = game->uap_vars;
        for (int i=0; i<game->uap_count; i++) {
            uint32_t bddvar = mtbdd_set_first(_vars);
            _vars = mtbdd_set_next(_vars);
            var_to_lit.emplace(bddvar, uap_to_lit[i]);
        }
    }

    // compute the set of states so we can create state_to_lit
    MTBDD states = sylvan_project(game->trans, game->s_vars);
    mtbdd_protect(&states);

    // make state_to_lit
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
            state_to_lit[state] = lit;
            lit += 2;
            // next state
            lf = mtbdd_enum_all_next(states, game->s_vars, state_arr, NULL);
        }
    }

    // data structures
    std::vector<std::vector<int>> uap_products;
    std::vector<std::vector<int>> uap_sums;
    std::vector<std::vector<int>> state_sums;
    std::vector<std::vector<int>> uap_state_products;
    std::vector<std::vector<int>> uap_state_sums;

    // do each controllable AP (output signals)
    MTBDD cap_bdd = mtbdd_false;
    mtbdd_protect(&cap_bdd);

    MTBDD s = mtbdd_false;
    mtbdd_protect(&s);

    for (int i=0; i<game->cap_count; i++) {
        cap_bdd = sylvan_ithvar(mtbdd_set_first(game->cap_vars)+i);
        // keep just s and u... get rid of other cap variables
        cap_bdd = sylvan_and_exists(game->strategies, cap_bdd, game->cap_vars);

        // this gives source states and UAP for this cap
        auto uaps = BDDTools::collectSubroots(cap_bdd, mtbdd_set_first(game->uns_vars));

        uap_state_sums.emplace_back();

        for (MTBDD uap : uaps) {
            s = BDDTools::pathsToSubroot(cap_bdd, mtbdd_set_first(game->uns_vars), uap);

            state_sums.emplace_back();
            uint8_t state_arr_2[game->statebits];
            MTBDD lf2 = mtbdd_enum_all_first(s, game->s_vars, state_arr_2, NULL);
            while (lf2 != mtbdd_false) {
                // decode state
                int from = 0;
                for (int i=0; i<game->statebits; i++) {
                    from <<= 1;
                    if (state_arr_2[i]) from |= 1;
                }
                // if state 0, take inverse
                state_sums.back().push_back(from == 0 ? aiger_not(state_to_lit.at(from)) : state_to_lit.at(from));
                lf2 = mtbdd_enum_all_next(s, game->s_vars, state_arr_2, NULL);
            }

            uap_sums.emplace_back();
            MTBDD bddres;
            ZDD isop = zdd_isop(uap, uap, &bddres);
            zdd_protect(&isop);
            // assert(uap == bddres);
            // assert(uap == zdd_cover_to_bdd(isop));
            if (isop == zdd_true) {
                // it is aiger_true
                uap_sums.back().push_back(-1);
            } else if (isop == zdd_false) {
                // it is aiger_false
                uap_sums.back().push_back(-2);
            } else {
                uap_idx = uap_products.size();
                // loop over all products
                int product[game->statebits+game->uap_count+1] = { 0 };
                ZDD res = zdd_cover_enum_first(isop, product);
                while (res != zdd_false) {
                    uap_products.emplace_back();
                    for (int i=0; product[i] != -1; i++) {
                        int the_lit = var_to_lit.at(product[i]/2);
                        if (product[i]&1) the_lit = aiger_not(the_lit);
                        uap_products.back().push_back(the_lit);
                    }
                    res = zdd_cover_enum_next(isop, product);
                    uap_sums.back().push_back(uap_products.size()-1);
                }
            }
            zdd_unprotect(&isop);

            uap_state_products.emplace_back();
            uap_state_products.back().push_back(uap_sums.size()-1);
            uap_state_products.back().push_back(state_sums.size()-1);
            uap_state_sums.back().push_back(uap_state_products.size()-1);
        }
    }

    MTBDD su_vars = mtbdd_set_addall(game->p_vars, game->cap_vars);
    mtbdd_protect(&su_vars);

    // full is: s > u > ns
    MTBDD full = mtbdd_and_exists(game->strategies, game->trans, su_vars);
    mtbdd_protect(&full);

    std::vector<int> states_vec;

    uint8_t state_arr[game->statebits];
    MTBDD lf = mtbdd_enum_all_first(states, game->s_vars, state_arr, NULL);
    while (lf != mtbdd_false) {
        // decode state
        int state = 0;
        for (int i=0; i<game->statebits; i++) {
            state <<= 1;
            if (state_arr[i]) state |= 1;
        }
        states_vec.push_back(state);

        // encode state using NS vars
        cap_bdd = encode_state(state, game->statebits, mtbdd_set_first(game->ns_vars));
        // keep s > u of this state
        cap_bdd = sylvan_and_exists(full, cap_bdd, game->ns_vars);

        // this gives source states and UAP for this state
        auto uaps = BDDTools::collectSubroots(cap_bdd, mtbdd_set_first(game->uns_vars));

        uap_state_sums.emplace_back();
        for (MTBDD uap : uaps) {
            s = BDDTools::pathsToSubroot(cap_bdd, mtbdd_set_first(game->uns_vars), uap);

            state_sums.emplace_back();
            uint8_t state_arr_2[game->statebits];
            MTBDD lf2 = mtbdd_enum_all_first(s, game->s_vars, state_arr_2, NULL);
            while (lf2 != mtbdd_false) {
                // decode state
                int from = 0;
                for (int i=0; i<game->statebits; i++) {
                    from <<= 1;
                    if (state_arr_2[i]) from |= 1;
                }
                // if state 0, take inverse
                state_sums.back().push_back(from == 0 ? aiger_not(state_to_lit.at(from)) : state_to_lit.at(from));
                lf2 = mtbdd_enum_all_next(s, game->s_vars, state_arr_2, NULL);
            }

            uap_sums.emplace_back();
            MTBDD bddres;
            ZDD isop = zdd_isop(uap, uap, &bddres);
            zdd_protect(&isop);
            // assert(uap == bddres);
            // assert(uap == zdd_cover_to_bdd(isop));
            if (isop == zdd_true) {
                // it is aiger_true
                uap_sums.back().push_back(-1);
            } else if (isop == zdd_false) {
                // it is aiger_false
                uap_sums.back().push_back(-2);
            } else {
                uap_idx = uap_products.size();
                // loop over all products
                int product[game->statebits+game->uap_count+1] = { 0 };
                ZDD res = zdd_cover_enum_first(isop, product);
                while (res != zdd_false) {
                    uap_products.emplace_back();
                    for (int i=0; product[i] != -1; i++) {
                        int the_lit = var_to_lit.at(product[i]/2);
                        if (product[i]&1) the_lit = aiger_not(the_lit);
                        uap_products.back().push_back(the_lit);
                    }
                    uap_sums.back().push_back(uap_products.size()-1);
                    res = zdd_cover_enum_next(isop, product);
                }
            }
            zdd_unprotect(&isop);

            uap_state_products.emplace_back();
            uap_state_products.back().push_back(uap_sums.size()-1);
            uap_state_products.back().push_back(state_sums.size()-1);
            uap_state_sums.back().push_back(uap_state_products.size()-1);
        }
        lf = mtbdd_enum_all_next(states, game->s_vars, state_arr, NULL);
    }

    mtbdd_unprotect(&states);
    mtbdd_unprotect(&cap_bdd);
    mtbdd_unprotect(&s);
    mtbdd_unprotect(&full);
    mtbdd_unprotect(&su_vars);

    // Process UAP_PRODUCTS
    if (verbose) std::cerr << "encoding input products" << std::endl << std::flush;
    reduce(uap_products, false);

    // Process UAP_SUMS
    if (verbose) std::cerr << "encoding input sums of products" << std::endl << std::flush;
    for (auto& x : uap_sums) {
        for (int i=0; i<x.size(); i++) {
            if (x[i] == -2) x[i] = aiger_false;
            else if (x[i] == -1) x[i] = aiger_true;
            else x[i] = uap_products[x[i]].front();
        }
    }
    reduce(uap_sums, true);

    // Process STATES_SUMS
    if (verbose) std::cerr << "encoding state sums" << std::endl << std::flush;
    reduce(state_sums, true);

    // Process UAP_STATE_PRODUCTS
    if (verbose) std::cerr << "encoding uap-state products" << std::endl << std::flush;
    for (auto& x : uap_state_products) {
        auto u = uap_sums[x[0]].front();
        auto s = state_sums[x[1]].front();
        x[0] = u;
        x[1] = s;
        // actually, skip reduce
        // we immediately make AND(u, s) instead of letting reduce do its magic
        // the reason being that this is very unlikely to find improvements anyhow and it takes "long"
        // (4 seconds on full_arbiter_8 which isn't so bad, but it did not improve the number of gates
        x.clear();
        x.push_back(makeand(u, s));
    }
    // reduce(uap_state_products, false);

    // Process UAP_STATE_SUMS
    if (verbose) std::cerr << "encoding uap-state sums of products" << std::endl << std::flush;
    for (auto& x : uap_state_sums) {
        if (x.size() == 0) {
            x.push_back(aiger_false);
        } else {
            for (int i=0; i<x.size(); i++) {
                x[i] = uap_state_products[x[i]].front();
            }
        }
    }
    reduce(uap_state_sums, true);

    for (int i=0; i<game->cap_count; i++) {
        // get the aiger thing
        auto result = uap_state_sums[i].front();
        aiger_add_output(a, result, caps[i]); // simple, really
    }

    for (int i=0; i<states_vec.size(); i++) {
        auto state = states_vec[i];
        auto result = uap_state_sums[game->cap_count + i].front();
        if (state == 0) result = aiger_not(result);
        aiger_add_latch(a, state_to_lit.at(state), result, "");
    }
}

void
AIGmaker::process()
{
    if (sop) return process_sop();

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
            var_to_lit.emplace(bddvar, uap_to_lit[i]);
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

            // this gives source states and UAP for this cap
            auto uaps = BDDTools::collectSubroots(cap_bdd, mtbdd_set_first(game->uns_vars));
            // each uap is a MTBDD node in cap_bdd, so no need to reference

            std::deque<int> terms;

            for (MTBDD uap : uaps) {
                // s is all states that go to that particular uap
                s = BDDTools::pathsToSubroot(cap_bdd, mtbdd_set_first(game->uns_vars), uap);

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
                    source_gates.push_back(from == 0 ? aiger_not(state_to_lit.at(from)) : state_to_lit.at(from));
                    lf2 = mtbdd_enum_all_next(s, game->s_vars, state_arr_2, NULL);
                }

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

            auto uaps = BDDTools::collectSubroots(cap_bdd, mtbdd_set_first(game->uns_vars));

            std::deque<int> terms;
            for (MTBDD uap : uaps) {
                s = BDDTools::pathsToSubroot(cap_bdd, mtbdd_set_first(game->uns_vars), uap);

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
                    source_gates.push_back(from == 0 ? aiger_not(state_to_lit.at(from)) : state_to_lit.at(from));
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
            // if state 0, take inverse
            if (state == 0) result = aiger_not(result);
            aiger_add_latch(a, state_to_lit.at(state), result, "");
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
            var_to_lit.emplace(bddvar, state_to_lit.at(i));
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
            }

            ZDD* state_zdds = new ZDD[game->statebits];
            for (int i=0; i<game->statebits; i++) {
                state_zdds[i] = zdd_false;
                zdd_protect(&state_zdds[i]);

                MTBDD bddres;
                state_zdds[i] = zdd_isop(state_bdds[i], state_bdds[i], &bddres);
                assert(bddres == state_bdds[i]);
            }

            for (int i=0; i<game->cap_count; i++) {
                int res = bdd_to_aig_cover_sop(cap_zdds[i]);
                aiger_add_output(a, res, caps[i]); // simple, really
            }
            for (int i=0; i<game->statebits; i++) {
                int res = bdd_to_aig_cover_sop(state_zdds[i]);
                aiger_add_latch(a, state_to_lit.at(i), res, "");
            }
        } else {
            for (int i=0; i<game->cap_count; i++) {
                int res = bdd_to_aig(cap_bdds[i]);
                aiger_add_output(a, res, caps[i]); // simple, really
            }
            for (int i=0; i<game->statebits; i++) {
                int res = bdd_to_aig(state_bdds[i]);
                aiger_add_latch(a, state_to_lit.at(i), res, "");
            }
        }

        for (int i=0; i<game->cap_count; i++) mtbdd_unprotect(&cap_bdds[i]);
        for (int i=0; i<game->statebits; i++) mtbdd_unprotect(&state_bdds[i]);
        delete[] cap_bdds;
        delete[] state_bdds;
    }
}

int
AIGmaker::writeAscii(FILE* out)
{
    return aiger_write_to_file(a, aiger_ascii_mode, out);
}

int
AIGmaker::writeBinary(FILE* out)
{
    return aiger_write_to_file(a, aiger_binary_mode, out);
}

void AIGmaker::readFile(FILE* infile) {
    aiger_reset(a);
    a = aiger_init();
    const char *read_result = aiger_read_from_file(a, infile);
    if (read_result != nullptr) {
        throw std::runtime_error("Could not read AIGER circuit from file: " + std::string(read_result));
    }
    aiger_delete_comments(a);
}

void AIGmaker::drewrite() {
    ABCMinimization min(*this, verbose);
    min.drewrite();
}

void AIGmaker::compress() {
    ABCMinimization min(*this, verbose);
    min.compress();
}