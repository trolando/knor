/**************************************************************************
 * Copyright Tom van Dijk
 *************************************************************************/

#include <abcminimization.hpp>
#include <aigencoder.hpp>
#include <bddtools.hpp>
#include <knor.hpp>
#include <algorithm>
#include <set>

using namespace sylvan;

AIGEncoder::AIGEncoder(HoaData &data, SymGame &game) : data(data), game(game), circuit(nullptr)
{
}

unsigned int AIGEncoder::bddToAigIsop(MTBDD bdd)
{
    if (verbose) {
        // std::cerr << "running isop for BDD with " << mtbdd_nodecount(bdd) << " nodes." << std::endl;
    }
    MTBDD bddres;
    ZDD the_isop = zdd_isop(bdd, bdd, &bddres);
    zdd_protect(&the_isop);
    // no need to reference the result...
    assert(bdd == bddres);
    assert(bdd == zdd_cover_to_bdd(the_isop));
    if (verbose) {
        // std::cerr << "isop has " << (long)zdd_pathcount(isop) << " terms and " << zdd_nodecount(&isop, 1) << " nodes." << std::endl;
    }

    auto res = bddToAigCoverSop(the_isop);
    zdd_unprotect(&the_isop);
    return res;
}

unsigned int AIGEncoder::bddToAigCoverSop(ZDD cover)
{
    if (cover == zdd_true) return aiger_true;
    if (cover == zdd_false) return aiger_false;

    // a queue that stores all products, which will need to be summed
    std::deque<unsigned int> products;

    // a product could consist of all variables, and a -1 to denote the end of the product
    int32_t product[game.statebits + game.uap_count + 1];
    auto res = zdd_cover_enum_first(cover, product);
    while (res != zdd_false) {
        //  containing subproducts in the form of gates
        std::deque<unsigned int> gates;

        for (int i=0; product[i] != -1; i++) {
            auto the_lit = bddvar_to_lit.at(product[i] / 2);
            if (product[i]&1) the_lit = aiger_not(the_lit);
            gates.push_back(the_lit);
        }

        // simplifyAnd(gates);

        // while we still have subproducts we need to AND together
        while (!gates.empty()) {
            auto last = gates.front();
            gates.pop_front();
            if (!gates.empty()) {
                auto last2 = gates.front();
                gates.pop_front();
                auto new_gate = circuit->makeAnd(last, last2);
                gates.push_back(new_gate);
            } else {
                products.push_back(last);
            }
        }
        res = zdd_cover_enum_next(cover, product); // go to the next product
    }

    // products queue should now be full of complete products that need to be summed

    // simplifyOr(products);

    while (products.size() > 1) {
        auto product1 = products.front();
        products.pop_front();
        auto product2 = products.front();
        products.pop_front();
        auto summed_product = aiger_not(circuit->makeAnd(aiger_not(product1), aiger_not(product2)));
        products.push_back(summed_product);
    }

    return products.front();
}

unsigned int AIGEncoder::bddToAigCover(ZDD cover)
{
    // TODO can we make this nonrecursive?
    if (cover == zdd_true) return aiger_true;
    if (cover == zdd_false) return aiger_false;

    auto the_var = zdd_getvar(cover);
    auto the_lit = bddvar_to_lit.at(the_var / 2);
    if (the_var & 1) the_lit = aiger_not(the_lit);

    auto low = zdd_getlow(cover);
    auto high = zdd_gethigh(cover);

    auto res = the_lit;

    if (high != zdd_true) {
        auto x = bddToAigCover(high);
        res = circuit->makeAnd(res, x);
    }

    if (low != zdd_false) {
        auto x = bddToAigCover(low);
        res = aiger_not(circuit->makeAnd(aiger_not(res), aiger_not(x)));
    }

    return res;
}

unsigned int AIGEncoder::bddToAigRecursive(MTBDD bdd)
{
    if (bdd == mtbdd_true) return aiger_true;
    if (bdd == mtbdd_false) return aiger_false;

    auto comp = false;
    if (bdd & sylvan_complement) {
        bdd ^= sylvan_complement;
        comp = true;
    }

    auto the_lit = bddvar_to_lit.at(mtbdd_getvar(bdd));

    MTBDD low = mtbdd_getlow(bdd);
    MTBDD high = mtbdd_gethigh(bdd);

    unsigned int res;

    if (low == mtbdd_false) {
        // only high (value 1)
        if (high == mtbdd_true) {
            // actually this is the end, just the lit
            res = the_lit;
        } else {
            // AND(the_lit, ...)
            auto rhs0 = the_lit;
            auto rhs1 = bddToAigRecursive(high);
            res = circuit->makeAnd(rhs0, rhs1);
        }
    } else if (high == mtbdd_false) {
        // only low (value 0)
        if (low == mtbdd_true) {
            // actually this is the end, just the lit, negated
            res = aiger_not(the_lit);
        } else {
            // AND(not the_lit, ...)
            auto rhs0 = aiger_not(the_lit);
            auto rhs1 = bddToAigRecursive(low);
            res = circuit->makeAnd(rhs0, rhs1);
        }
    } else {
        // OR(low, high) == ~AND(~AND(the_lit, ...), ~AND(~the_lit, ...))
        auto lowres = bddToAigRecursive(low);
        auto highres = bddToAigRecursive(high);
        auto rhs0 = aiger_not(circuit->makeAnd(aiger_not(the_lit), lowres));
        auto rhs1 = aiger_not(circuit->makeAnd(the_lit, highres));
        res = aiger_not(circuit->makeAnd(rhs0, rhs1));
    }

    return comp ? aiger_not(res) : res;
}

void AIGEncoder::reduce(std::vector<std::vector<unsigned int>>& system, bool is_or)
{
    std::map<unsigned int, int> pop;
    for (auto& x : system) {
        if (x.size() == 1) continue;
        for (auto& y : x) {
            pop[y] += 1;
        }
    }

    std::map<unsigned int, int> pop2;

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
        auto m = 0, n = 0;
        auto gate = is_or ? aiger_not(circuit->makeAnd(aiger_not(first), aiger_not(second))) : circuit->makeAnd(first,
                                                                                                                second);
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
AIGEncoder::processSOP()
{
    // compute the set of states, so we can create state_to_lit
    MTBDD states = sylvan_project(game.strategies, game.s_vars);
    mtbdd_protect(&states);

    // make state_to_lit
    {
        uint8_t state_arr[game.statebits];
        MTBDD lf = mtbdd_enum_all_first(states, game.s_vars, state_arr, nullptr);
        while (lf != mtbdd_false) {
            // decode state
            int state = 0;
            for (int i=0; i<game.statebits; i++) {
                state <<= 1;
                if (state_arr[i]) state |= 1;
            }
            // give the state a literal
            state_to_lit[state] = circuit->makeLatch();
            // next state
            lf = mtbdd_enum_all_next(states, game.s_vars, state_arr, nullptr);
        }
    }

    // data structures
    std::vector<std::vector<unsigned int>> uap_products;
    std::vector<std::vector<unsigned int>> uap_sums;
    std::vector<std::vector<unsigned int>> state_sums;
    std::vector<std::vector<unsigned int>> uap_state_products;
    std::vector<std::vector<unsigned int>> uap_state_sums;

    MTBDD cap_bdd = mtbdd_false;
    mtbdd_protect(&cap_bdd);

    MTBDD s = mtbdd_false;
    mtbdd_protect(&s);

    // do each controllable AP (output signals)

    for (int i=0; i<game.cap_count; i++) {
        // TODO: remove assumption that the CAP vars are consecutive
        cap_bdd = sylvan_ithvar(mtbdd_set_first(game.cap_vars) + i);
        // keep just s and u... get rid of other cap variables
        cap_bdd = sylvan_and_exists(game.strategies, cap_bdd, game.cap_vars);

        {
            uint8_t state_arr_2[game.statebits];
            MTBDD lf2 = mtbdd_enum_all_first(cap_bdd, game.s_vars, state_arr_2, nullptr);
            while (lf2 != mtbdd_false) {
                // decode state
                int from = 0;
                for (int k=0; k<game.statebits; k++) {
                    from <<= 1;
                    if (state_arr_2[k]) from |= 1;
                }
                lf2 = mtbdd_enum_all_next(cap_bdd, game.s_vars, state_arr_2, nullptr);
            }
        }

        // this gives source states and UAP for this cap
        auto uaps = BDDTools::collectSubroots(cap_bdd, mtbdd_set_first(game.uns_vars));

        uap_state_sums.emplace_back();

        for (MTBDD uap : uaps) {
            s = BDDTools::pathsToSubroot(cap_bdd, mtbdd_set_first(game.uns_vars), uap);

            state_sums.emplace_back();

            uint8_t state_arr_2[game.statebits];
            MTBDD lf2 = mtbdd_enum_all_first(s, game.s_vars, state_arr_2, nullptr);
            while (lf2 != mtbdd_false) {
                // decode state
                int from = 0;
                for (int k=0; k<game.statebits; k++) {
                    from <<= 1;
                    if (state_arr_2[k]) from |= 1;
                }
                // if state 0, take inverse
                state_sums.back().push_back(from == 0 ? aiger_not(state_to_lit.at(from)) : state_to_lit.at(from));
                lf2 = mtbdd_enum_all_next(s, game.s_vars, state_arr_2, nullptr);
            }

            uap_sums.emplace_back();
            MTBDD bddres;
            ZDD the_isop = zdd_isop(uap, uap, &bddres);
            zdd_protect(&the_isop);
            // assert(uap == bddres);
            // assert(uap == zdd_cover_to_bdd(the_isop));
            if (the_isop == zdd_true) {
                // it is aiger_true
                uap_sums.back().push_back(-1);
            } else if (the_isop == zdd_false) {
                // it is aiger_false
                uap_sums.back().push_back(-2);
            } else {
                // loop over all products
                int product[game.statebits + game.uap_count + 1];
                ZDD res = zdd_cover_enum_first(the_isop, product);
                while (res != zdd_false) {
                    uap_products.emplace_back();
                    for (int k=0; product[k] != -1; k++) {
                        auto the_lit = bddvar_to_lit.at(product[k] / 2);
                        if (product[k]&1) the_lit = aiger_not(the_lit);
                        uap_products.back().push_back(the_lit);
                    }
                    res = zdd_cover_enum_next(the_isop, product);
                    uap_sums.back().push_back(uap_products.size()-1);
                }
            }
            zdd_unprotect(&the_isop);

            uap_state_products.emplace_back();
            uap_state_products.back().push_back(uap_sums.size()-1);
            uap_state_products.back().push_back(state_sums.size()-1);
            uap_state_sums.back().push_back(uap_state_products.size()-1);
        }
    }

    MTBDD su_vars = mtbdd_set_addall(game.np_vars, game.cap_vars);
    mtbdd_protect(&su_vars);

    // full is: s > u > ns
    MTBDD full = mtbdd_and_exists(game.strategies, game.trans, su_vars);
    mtbdd_protect(&full);

    std::vector<int> states_vec;

    uint8_t state_arr[game.statebits];
    MTBDD lf = mtbdd_enum_all_first(states, game.s_vars, state_arr, nullptr);
    while (lf != mtbdd_false) {
        // decode state
        auto state = 0;
        for (int i=0; i<game.statebits; i++) {
            state <<= 1;
            if (state_arr[i]) state |= 1;
        }
        states_vec.push_back(state);

        // encode state using NS vars
        cap_bdd = BDDTools::encode_state(state, game.ns_vars);
        // keep s > u of this state
        cap_bdd = sylvan_and_exists(full, cap_bdd, game.ns_vars);

        // this gives source states and UAP for this state
        auto uaps = BDDTools::collectSubroots(cap_bdd, mtbdd_set_first(game.uns_vars));

        uap_state_sums.emplace_back();
        for (MTBDD uap : uaps) {
            s = BDDTools::pathsToSubroot(cap_bdd, mtbdd_set_first(game.uns_vars), uap);

            state_sums.emplace_back();
            uint8_t state_arr_2[game.statebits];
            MTBDD lf2 = mtbdd_enum_all_first(s, game.s_vars, state_arr_2, nullptr);
            while (lf2 != mtbdd_false) {
                // decode state
                auto from = 0;
                for (int i=0; i<game.statebits; i++) {
                    from <<= 1;
                    if (state_arr_2[i]) from |= 1;
                }
                // if state 0, take inverse
                state_sums.back().push_back(from == 0 ? aiger_not(state_to_lit.at(from)) : state_to_lit.at(from));
                lf2 = mtbdd_enum_all_next(s, game.s_vars, state_arr_2, nullptr);
            }

            uap_sums.emplace_back();
            MTBDD bddres;
            auto the_isop = zdd_isop(uap, uap, &bddres);
            zdd_protect(&the_isop);
            // assert(uap == bddres);
            // assert(uap == zdd_cover_to_bdd(the_isop));
            if (the_isop == zdd_true) {
                // it is aiger_true
                uap_sums.back().push_back(-1);
            } else if (the_isop == zdd_false) {
                // it is aiger_false
                uap_sums.back().push_back(-2);
            } else {
                // loop over all products
                auto product = std::vector<int>(game.statebits + game.uap_count + 1, 0);
                auto res = zdd_cover_enum_first(the_isop, product.data());
                while (res != zdd_false) {
                    uap_products.emplace_back();
                    for (int i=0; product[i] != -1; i++) {
                        auto the_lit = bddvar_to_lit.at(product[i] / 2);
                        if (product[i]&1) the_lit = aiger_not(the_lit);
                        uap_products.back().push_back(the_lit);
                    }
                    uap_sums.back().push_back(uap_products.size() - 1);
                    res = zdd_cover_enum_next(the_isop, product.data());
                }
            }
            zdd_unprotect(&the_isop);

            uap_state_products.emplace_back();
            uap_state_products.back().push_back(uap_sums.size()-1);
            uap_state_products.back().push_back(state_sums.size()-1);
            uap_state_sums.back().push_back(uap_state_products.size()-1);
        }
        lf = mtbdd_enum_all_next(states, game.s_vars, state_arr, nullptr);
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
        for (unsigned int & i : x) {
            if (i == -2) i = aiger_false;
            else if (i == -1) i = aiger_true;
            else i = uap_products[i].front();
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
        auto sum = state_sums[x[1]].front();
        x[0] = u;
        x[1] = sum;
        // actually, skip reduce
        // we immediately make AND(u, s) instead of letting reduce do its magic
        // the reason being that this is very unlikely to find improvements anyhow and it takes "long"
        // (4 seconds on full_arbiter_8 which isn't so bad, but it did not improve the number of gates
        x.clear();
        x.push_back(circuit->makeAnd(u, sum));
    }
    // reduce(uap_state_products, false);

    // Process UAP_STATE_SUMS
    if (verbose) std::cerr << "encoding uap-state sums of products" << std::endl << std::flush;
    for (auto& x : uap_state_sums) {
        if (x.empty()) {
            x.push_back(aiger_false);
        } else {
            for (unsigned int & i : x) {
                i = uap_state_products[i].front();
            }
        }
    }
    reduce(uap_state_sums, true);

    for (int i=0; i<game.cap_count; i++) {
        // get the aiger thing
        auto result = uap_state_sums[i].front();
        circuit->makeOutput(result, cap_labels[i]);
    }

    for (int i=0; i<states_vec.size(); i++) {
        auto state = states_vec[i];
        auto result = uap_state_sums[game.cap_count + i].front();
        if (state == 0) result = aiger_not(result);
        circuit->setLatch(state_to_lit.at(state), result, "");
    }
}

void AIGEncoder::processOnehot()
{
    MTBDD states = sylvan_project(game.trans, game.s_vars);
    mtbdd_protect(&states);

    /* prepare state_to_lit */
    {
        uint8_t state_arr[game.statebits];
        MTBDD lf = mtbdd_enum_all_first(states, game.s_vars, state_arr, nullptr);
        while (lf != mtbdd_false) {
            // decode state
            int state = 0;
            for (int i=0; i<game.statebits; i++) {
                state <<= 1;
                if (state_arr[i]) state |= 1;
            }
            // give the state a literal
            // std::cerr << "state " << state << " gets latch literal " << lit << std::endl;
            state_to_lit[state] = circuit->makeLatch();
            // next state
            lf = mtbdd_enum_all_next(states, game.s_vars, state_arr, nullptr);
        }
    }

    // do each controllable AP (output signals)
    MTBDD cap_bdd = mtbdd_false;
    mtbdd_protect(&cap_bdd);

    MTBDD s = mtbdd_false;
    mtbdd_protect(&s);

    for (int i=0; i<game.cap_count; i++) {
        cap_bdd = sylvan_ithvar(mtbdd_set_first(game.cap_vars)+i);
        // keep just s and u... get rid of other cap variables
        cap_bdd = sylvan_and_exists(game.strategies, cap_bdd, game.cap_vars);

        // this gives source states and UAP for this cap
        auto uaps = BDDTools::collectSubroots(cap_bdd, mtbdd_set_first(game.uns_vars));
        // each uap is a MTBDD node in cap_bdd, so no need to reference

        std::deque<unsigned int> terms;

        for (MTBDD uap : uaps) {
            // s is all states that go to that particular uap
            s = BDDTools::pathsToSubroot(cap_bdd, mtbdd_set_first(game.uns_vars), uap);

            std::vector<int> source_states;
            std::deque<unsigned int> source_gates;

            uint8_t state_arr_2[game.statebits];
            MTBDD lf2 = mtbdd_enum_all_first(s, game.s_vars, state_arr_2, nullptr);
            while (lf2 != mtbdd_false) {
                // decode state
                auto from = 0;
                for (auto k=0; k<game.statebits; k++) {
                    from <<= 1;
                    if (state_arr_2[k]) from |= 1;
                }
                source_states.push_back(from);
                // if state 0, take inverse
                source_gates.push_back(from == 0 ? aiger_not(state_to_lit.at(from)) : state_to_lit.at(from));
                lf2 = mtbdd_enum_all_next(s, game.s_vars, state_arr_2, nullptr);
            }

            auto aig_uap = isop ? bddToAigIsop(uap) : bddToAigRecursive(uap);
            auto aig_states = circuit->makeMultiOr(source_gates);
            terms.push_back(circuit->makeAnd(aig_uap, aig_states));
        }

        auto result = circuit->makeMultiOr(terms);
        circuit->makeOutput(result, cap_labels[i]);
    }

    MTBDD su_vars = mtbdd_set_addall(game.np_vars, game.cap_vars);
    mtbdd_protect(&su_vars);

    // full is: s > u > ns
    MTBDD full = mtbdd_and_exists(game.strategies, game.trans, su_vars);
    mtbdd_protect(&full);

    // long noStates = (long)sylvan_satcount(states, game.s_vars);
    // std::cerr << "number of states: " << noStates << std::endl;

    uint8_t state_arr[game.statebits];
    MTBDD lf = mtbdd_enum_all_first(states, game.s_vars, state_arr, nullptr);
    while (lf != mtbdd_false) {
        // decode state
        int state = 0;
        for (int i=0; i<game.statebits; i++) {
            state <<= 1;
            if (state_arr[i]) state |= 1;
        }

        // encode state using NS vars
        cap_bdd = BDDTools::encode_state(state, game.ns_vars);
        // keep s > u of this state
        cap_bdd = sylvan_and_exists(full, cap_bdd, game.ns_vars);

        // std::cerr << "state " << state << /*" to has " << s <<*/ std::endl;

        auto uaps = BDDTools::collectSubroots(cap_bdd, mtbdd_set_first(game.uns_vars));

        std::deque<unsigned int> terms;
        for (MTBDD uap : uaps) {
            s = BDDTools::pathsToSubroot(cap_bdd, mtbdd_set_first(game.uns_vars), uap);

            std::vector<int> source_states;
            std::deque<unsigned int> source_gates;

            uint8_t state_arr_2[game.statebits];
            MTBDD lf2 = mtbdd_enum_all_first(s, game.s_vars, state_arr_2, nullptr);
            while (lf2 != mtbdd_false) {
                // decode state
                int from = 0;
                for (int i=0; i<game.statebits; i++) {
                    from <<= 1;
                    if (state_arr_2[i]) from |= 1;
                }
                source_states.push_back(from);
                // if state 0, take inverse
                source_gates.push_back(from == 0 ? aiger_not(state_to_lit.at(from)) : state_to_lit.at(from));
                lf2 = mtbdd_enum_all_next(s, game.s_vars, state_arr_2, nullptr);
            }

            //for (int from : source_states) {
            // std::cerr << "source state " << from << " with UAP " << ((uap&sylvan_complement)?"~":"") << (uap&~sylvan_complement) << std::endl;
            //}

            auto aig_uap = isop ? bddToAigIsop(uap) : bddToAigRecursive(uap);
            auto aig_states = circuit->makeMultiOr(source_gates);
            auto result = circuit->makeAnd(aig_uap, aig_states);
            terms.push_back(result);
        }

        auto result = circuit->makeMultiOr(terms);
        // if state 0, take inverse
        if (state == 0) result = aiger_not(result);
        circuit->setLatch(state_to_lit.at(state), result, "");
        lf = mtbdd_enum_all_next(states, game.s_vars, state_arr, nullptr);
    }

    mtbdd_unprotect(&states);
    mtbdd_unprotect(&cap_bdd);
    mtbdd_unprotect(&s);
    mtbdd_unprotect(&full);
    mtbdd_unprotect(&su_vars);
}

void AIGEncoder::processBinary()
{
    // initialize state_to_lit and bddvar_to_lit
    {
        auto _vars = game.s_vars;
        for (int i = 0; i < game.statebits; i++) {
            auto bddvar = mtbdd_set_first(_vars);
            _vars = mtbdd_set_next(_vars);
            auto lit = circuit->makeLatch();
            state_to_lit[i] = lit;
            bddvar_to_lit.emplace(bddvar, lit);
        }
    }

    MTBDD cap_bdds[game.cap_count];   // contains the solution: controllable ap bdds: state -> uap -> B
    MTBDD state_bdds[game.statebits]; // contains the solution: state bit bdds      : state -> uap -> B

    // we use the "strategies" field of the game class, which should only include state -> uap -> cap.

    // compute the BDDs for the controllable APs
    for (auto i = 0; i < game.cap_count; i++) {
        cap_bdds[i] = sylvan_ithvar(mtbdd_set_first(game.cap_vars) + i);
        mtbdd_protect(&cap_bdds[i]);
        // keep just s and u... get rid of cap variables
        cap_bdds[i] = sylvan_and_exists(game.strategies, cap_bdds[i], game.cap_vars);
    }

    // now we compute the full BDD excluding priority and controllable AP
    // the BDD "full" is defined on variables state -> uap -> next state
    MTBDD pc_vars = mtbdd_set_addall(game.np_vars, game.cap_vars);
    mtbdd_protect(&pc_vars);
    MTBDD full = mtbdd_and_exists(game.strategies, game.trans, pc_vars);
    mtbdd_protect(&full);

    // compute bdds for each state variable
    for (int i=0; i<game.statebits; i++) {
        state_bdds[i] = sylvan_ithvar(mtbdd_set_first(game.ns_vars) + i);
        mtbdd_protect(&state_bdds[i]);
        // keep just s and u... get rid of next state variables
        state_bdds[i] = sylvan_and_exists(full, state_bdds[i], game.ns_vars);
    }

    mtbdd_unprotect(&full);
    mtbdd_unprotect(&pc_vars);

    if (isop) {
        // ISOP construction, first convert all cap bdds etc to covers
        ZDD cap_zdds[game.cap_count];
        for (int i=0; i<game.cap_count; i++) {
            cap_zdds[i] = zdd_false;
            zdd_protect(&cap_zdds[i]);

            MTBDD bddres;
            cap_zdds[i] = zdd_isop(cap_bdds[i], cap_bdds[i], &bddres);
            assert(bddres == cap_bdds[i]);
        }
        ZDD state_zdds[game.statebits];
        for (int i=0; i<game.statebits; i++) {
            state_zdds[i] = zdd_false;
            zdd_protect(&state_zdds[i]);

            MTBDD bddres;
            state_zdds[i] = zdd_isop(state_bdds[i], state_bdds[i], &bddres);
            assert(bddres == state_bdds[i]);
        }
        for (int i=0; i<game.cap_count; i++) {
            auto res = bddToAigCoverSop(cap_zdds[i]);
            circuit->makeOutput(res, cap_labels[i]);
            zdd_unprotect(&cap_zdds[i]);
        }
        for (int i=0; i<game.statebits; i++) {
            auto res = bddToAigCoverSop(state_zdds[i]);
            circuit->setLatch(state_to_lit.at(i), res, "");
            zdd_unprotect(&state_zdds[i]);
        }
    } else {
        // simple ITE (shannon expansion) construction
        if (false) {
            // EXPERIMENTAL ... reorder!
            // convert bdd vars to vector
            std::vector<unsigned int> bddvars;
            MTBDD _vars = game.s_vars;
            for (int i=0; i<game.statebits; i++) {
                bddvars.push_back(mtbdd_set_first(_vars));
                _vars = mtbdd_set_next(_vars);
            }
            // make the composition
            MTBDD s_to_s = mtbdd_map_empty();
            mtbdd_protect(&s_to_s);
            for (int i=0; i<game.statebits; i++) {
                s_to_s = mtbdd_map_add(s_to_s, bddvars.at(i), sylvan_ithvar(bddvars.at(game.statebits-i-1)));
            }
            for (int i=0; i<game.cap_count; i++) {
                cap_bdds[i] = sylvan_compose(cap_bdds[i], s_to_s);
            }
            for (int i=0; i<game.statebits; i++) {
                state_bdds[i] = sylvan_compose(state_bdds[i], s_to_s);
            }
            mtbdd_unprotect(&s_to_s);
            for (int i=0; i<game.statebits; i++) {
                bddvar_to_lit.insert_or_assign(bddvars.at(i), state_to_lit.at(game.statebits - i - 1));
            }
        }
        for (int i=0; i<game.cap_count; i++) {
            auto res = bddToAigRecursive(cap_bdds[i]);
            circuit->makeOutput(res, cap_labels[i]);
        }
        for (int i=0; i<game.statebits; i++) {
            auto res = bddToAigRecursive(state_bdds[i]);
            circuit->setLatch(state_to_lit.at(i), res, "");
        }
    }

    for (int i=0; i<game.cap_count; i++) mtbdd_unprotect(&cap_bdds[i]);
    for (int i=0; i<game.statebits; i++) mtbdd_unprotect(&state_bdds[i]);
}

std::unique_ptr<AIGCircuit> AIGEncoder::encode()
{
    circuit = std::make_unique<AIGCircuit>();

    // set which APs are controllable in the bitset controllable
    pg::bitset controllable(data.noAPs);
    for (int i=0; i<data.noCntAPs; i++) {
        controllable[data.cntAPs[i]] = true;
    }

    // fill the two arrays
    for (int i=0; i<data.noAPs; i++) {
        if (!controllable[i]) {
            uap_to_lit.push_back(circuit->makeInput(std::string(data.aps[i])));
        } else {
            cap_labels.emplace_back(data.aps[i]);
        }
    }

    // make bddvar_to_lit for uncontrollable APs
    MTBDD _vars = game.uap_vars;
    int index = 0;
    while (_vars != mtbdd_set_empty()) {
        uint32_t bddvar = mtbdd_set_first(_vars);
        _vars = mtbdd_set_next(_vars);
        bddvar_to_lit.emplace(bddvar, uap_to_lit[index++]);
    }

    if (sop) {
        processSOP();
    } else if (onehot) {
        processOnehot();
    } else {
        processBinary();
    }

    return std::move(circuit);
}
