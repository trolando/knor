/**************************************************************************
 * Copyright Tom van Dijk
 *************************************************************************/

#include <abcminimization.hpp>
#include <aigencoder.hpp>
#include <bddtools.hpp>
#include <knor.hpp>
#include <algorithm>
#include <set>
#include <unordered_map>

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
    } else if (high == mtbdd_true && rules.orGate) {
        // OR(the_lit, ...) == ~AND(~the_lit, ~...): one gate (the low edge is never complemented,
        // so it is never true)
        auto lowres = bddToAigRecursive(low);
        res = aiger_not(circuit->makeAnd(aiger_not(the_lit), aiger_not(lowres)));
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

// Choose the output and next-state functions one bit at a time among the
// allowed moves (state -> uap -> cap -> next state). Where a bit may be both 0
// and 1, or the state and input do not occur (such as unused state codes), its
// value is a don't care: any function between "1 where it must be 1" and "1
// unless it must be 0" will do. The smallest BDD of the former, its
// simplification by sylvan_restrict (on the care set) and the irredundant sum
// of products between the two is taken, with allCandidates also considering
// the latter and the simplification by sylvan_constrain. Each choice
// restricts the remaining moves, so together
// the functions select an allowed move in every state and for every input.
// The bits are chosen in the order outputs, state, or the reverse. The
// relation is the moves (state, inputVars -> outputs, next state); the
// functions are over its state and input variables.
void AIGEncoder::chooseWithDontCares(MTBDD relation, MTBDD inputVars, MTBDD* cap_bdds, MTBDD* state_bdds)
{
    MTBDD moves = mtbdd_false, out = mtbdd_false, bit = mtbdd_false, can1 = mtbdd_false, can0 = mtbdd_false;
    MTBDD fixed = mtbdd_false, must1 = mtbdd_false, candidate = mtbdd_false, link = mtbdd_false;
    for (MTBDD* r : {&moves, &out, &bit, &can1, &can0, &fixed, &must1, &candidate, &link}) mtbdd_refs_pushptr(r);
    moves = relation;
    out = mtbdd_set_addall(game.cap_vars, game.ns_vars);
    // Every state with a move must have one for every input: elsewhere the bits are don't cares, so a
    // state without a move for some input would get an arbitrary one.
    can1 = sylvan_exists(moves, out);
    can0 = sylvan_forall(can1, inputVars);
    fixed = sylvan_exists(can1, inputVars);
    if (can0 != fixed) {
        std::cerr << "error: a state of the controller has no move for some input" << std::endl;
        exit(1);
    }
    auto choose = [&](uint32_t var, MTBDD& result) {
        result = mtbdd_false;
        mtbdd_protect(&result);
        bit = sylvan_ithvar(var);
        can1 = sylvan_and_exists(moves, bit, out);
        can0 = sylvan_and_exists(moves, sylvan_not(bit), out);
        fixed = sylvan_xor(can1, can0);
        must1 = sylvan_and(can1, fixed);
        result = must1;
        auto consider = [&]() {
            if (mtbdd_nodecount(candidate) < mtbdd_nodecount(result)) result = candidate;
        };
        candidate = sylvan_restrict(must1, fixed);
        consider();
        candidate = sylvan_and(can0, fixed);        // must be 0
        candidate = sylvan_not(candidate);          // the upper bound
        if (allCandidates) consider();
        // the irredundant sum of products between the bounds (Minato-Morreale): its BDD is often
        // more regular than the generalized cofactors that restrict and constrain give
        if (isopCandidate) {
            zdd_isop(must1, candidate, &candidate);
            consider();
        }
        if (allCandidates) {
            candidate = sylvan_constrain(must1, fixed);
            consider();
        }
        link = sylvan_equiv(bit, result);
        moves = sylvan_and(moves, link);
    };
    auto outputs = [&]() {
        for (int i = 0; i < game.cap_count; i++) choose(mtbdd_set_first(game.cap_vars) + i, cap_bdds[i]);
    };
    auto states = [&]() {
        for (int i = 0; i < game.statebits; i++) choose(mtbdd_set_first(game.ns_vars) + i, state_bdds[i]);
    };
    if (statesFirst) {
        states();
        outputs();
    } else {
        outputs();
        states();
    }
    mtbdd_refs_popptr(9);
}

// The next-state function h, as ITE(g, h where g holds, h where g does not
// hold) for the output g whose logic the circuit already has, if the two
// restrictions (by sylvan_restrict, so this is exact) are smaller than h. For
// example a flag that a grant resets only needs the grant and a few gates.
unsigned int AIGEncoder::latchFunction(MTBDD h, MTBDD* cap_bdds, const std::vector<unsigned int>& output_lits)
{
    MTBDD high = mtbdd_false, low = mtbdd_false, bestHigh = mtbdd_false, bestLow = mtbdd_false;
    for (MTBDD* r : {&high, &low, &bestHigh, &bestLow}) mtbdd_refs_pushptr(r);
    size_t best = mtbdd_nodecount(h);
    int output = -1;
    for (int j = 0; j < game.cap_count; j++) {
        if (sylvan_isconst(cap_bdds[j])) continue;
        high = sylvan_restrict(h, cap_bdds[j]);
        low = sylvan_restrict(h, sylvan_not(cap_bdds[j]));
        const MTBDD both[2] = {high, low};
        const size_t cost = mtbdd_nodecount_more(both, 2) + 3;
        if (cost < best) {
            best = cost;
            output = j;
            bestHigh = high;
            bestLow = low;
        }
    }
    unsigned int res;
    if (output == -1) {
        res = bddToAigRecursive(h);
    } else {
        const unsigned int g = output_lits[output];
        const unsigned int a = bddToAigRecursive(bestHigh);
        const unsigned int b = bddToAigRecursive(bestLow);
        res = aiger_not(circuit->makeAnd(aiger_not(circuit->makeAnd(g, a)), aiger_not(circuit->makeAnd(aiger_not(g), b))));
    }
    mtbdd_refs_popptr(4);
    return res;
}

namespace {

// MTBDDs that stay protected as long as the pool exists; their addresses stay valid.
class Pool final {
public:
    Pool() = default;
    Pool(const Pool&) = delete;
    Pool& operator=(const Pool&) = delete;
    ~Pool() { for (MTBDD& dd : items) mtbdd_unprotect(&dd); }

    MTBDD& add(MTBDD dd)
    {
        items.push_back(dd);
        mtbdd_protect(&items.back());
        return items.back();
    }

private:
    std::deque<MTBDD> items;
};

// The functions of the inputs that the states use for one signal, and which states use each.
struct Functions {
    std::vector<MTBDD*> function;
    std::vector<std::vector<int>> states;
    std::map<MTBDD, int> index;
};

}  // namespace

// One-hot encoding that chooses among the allowed moves: a latch per state
// (that of the initial state inverted, so that all latches start at 0), and
// for every output and every latch an OR over the states of "in this state
// and this function of the inputs". Per state, the outputs and the successor
// are chosen one signal at a time, as in chooseWithDontCares but as functions
// of the inputs only, preferring a function that another state already uses
// for the same signal (which the circuit then shares). Every state chooses one
// successor for every input, so exactly one state is active at any time.
// Returns false (without a circuit) if the game has more than maxStates states.
bool AIGEncoder::processOnehotDontCares(int maxStates)
{
    Pool pool;
    MTBDD& moves = pool.add(sylvan_and_exists(game.strategies, game.trans, game.np_vars));   // s, u, c, ns
    MTBDD& out = pool.add(mtbdd_set_addall(game.cap_vars, game.ns_vars));
    MTBDD& states = pool.add(sylvan_project(moves, game.s_vars));
    if (sylvan_satcount(states, game.s_vars) > maxStates) return false;

    // the states, and per state its moves (u, c, ns)
    std::vector<int> code;
    std::map<int, int> index;   // state code -> position
    std::vector<MTBDD*> rel;
    {
        std::vector<uint8_t> arr(game.statebits);
        MTBDD& cube = pool.add(mtbdd_false);
        for (MTBDD leaf = mtbdd_enum_all_first(states, game.s_vars, arr.data(), nullptr); leaf != mtbdd_false;
             leaf = mtbdd_enum_all_next(states, game.s_vars, arr.data(), nullptr)) {
            int c = 0;
            for (int i = 0; i < game.statebits; i++) c = 2 * c + (arr[i] ? 1 : 0);
            index[c] = static_cast<int>(code.size());
            code.push_back(c);
        }
        for (int c : code) {
            cube = BDDTools::encode_state(c, game.s_vars);
            rel.push_back(&pool.add(sylvan_and_exists(moves, cube, game.s_vars)));
        }
    }
    if (index.count(0) == 0) return false;   // no initial state
    const int n = static_cast<int>(code.size());
    {
        // every successor must be a state with moves
        MTBDD& successors = pool.add(sylvan_project(moves, game.ns_vars));
        std::vector<uint8_t> arr(game.statebits);
        for (MTBDD leaf = mtbdd_enum_all_first(successors, game.ns_vars, arr.data(), nullptr); leaf != mtbdd_false;
             leaf = mtbdd_enum_all_next(successors, game.ns_vars, arr.data(), nullptr)) {
            int c = 0;
            for (int i = 0; i < game.statebits; i++) c = 2 * c + (arr[i] ? 1 : 0);
            if (index.count(c) == 0) return false;
        }
    }

    MTBDD& can1 = pool.add(mtbdd_false);
    MTBDD& can0 = pool.add(mtbdd_false);
    MTBDD& fixed = pool.add(mtbdd_false);
    MTBDD& must1 = pool.add(mtbdd_false);
    MTBDD& must0 = pool.add(mtbdd_false);
    MTBDD& candidate = pool.add(mtbdd_false);
    MTBDD& chosen = pool.add(mtbdd_false);
    MTBDD& test = pool.add(mtbdd_false);
    MTBDD& lit = pool.add(mtbdd_false);
    MTBDD& link = pool.add(mtbdd_false);
    MTBDD& rest = pool.add(mtbdd_false);

    // Choose a function of the inputs between must1 and not must0 for one signal of state k: a
    // constant, or a function that other states use for it (the most used first), or the smallest
    // of the simplifications; records it in f and returns it.
    auto choose = [&](Functions& f, int k) -> MTBDD {
        fixed = sylvan_xor(can1, can0);
        must1 = sylvan_and(can1, fixed);
        must0 = sylvan_and(can0, fixed);
        auto fits = [&](MTBDD g) {
            test = sylvan_and(must1, sylvan_not(g));
            if (test != mtbdd_false) return false;
            test = sylvan_and(g, must0);
            return test == mtbdd_false;
        };
        bool found = true;
        if (must1 == mtbdd_false) {
            chosen = mtbdd_false;
        } else if (must0 == mtbdd_false) {
            chosen = mtbdd_true;
        } else {
            found = false;
            std::vector<int> order(f.function.size());
            for (size_t j = 0; j < order.size(); j++) order[j] = static_cast<int>(j);
            std::sort(order.begin(), order.end(), [&](int a, int b) { return f.states[a].size() > f.states[b].size(); });
            for (size_t j = 0; j < order.size() && j < 16 && !found; j++) {
                if (fits(*f.function[order[j]])) {
                    chosen = *f.function[order[j]];
                    found = true;
                }
            }
            if (!found) {
                chosen = must1;
                auto consider = [&]() {
                    if (mtbdd_nodecount(candidate) < mtbdd_nodecount(chosen)) chosen = candidate;
                };
                candidate = sylvan_restrict(must1, fixed);
                consider();
                candidate = sylvan_not(must0);
                if (allCandidates) consider();
                if (isopCandidate) {
                    zdd_isop(must1, candidate, &candidate);
                    consider();
                }
                if (allCandidates) {
                    candidate = sylvan_constrain(must1, fixed);
                    consider();
                }
            }
        }
        if (chosen != mtbdd_false) {
            auto it = f.index.find(chosen);
            if (it == f.index.end()) {
                it = f.index.emplace(chosen, static_cast<int>(f.function.size())).first;
                f.function.push_back(&pool.add(chosen));
                f.states.emplace_back();
            }
            f.states[it->second].push_back(k);
        }
        return chosen;
    };

    std::vector<Functions> outputs(game.cap_count), latches(n);
    auto chooseOutputs = [&](int k) {
        MTBDD& r = *rel[k];
        for (int i = 0; i < game.cap_count; i++) {
            lit = sylvan_ithvar(mtbdd_set_first(game.cap_vars) + i);
            can1 = sylvan_and_exists(r, lit, out);
            can0 = sylvan_and_exists(r, sylvan_not(lit), out);
            const MTBDD f = choose(outputs[i], k);
            link = sylvan_equiv(lit, f);
            r = sylvan_and(r, link);
        }
    };
    auto chooseSuccessor = [&](int k) {
        MTBDD& r = *rel[k];
        // the successors, those allowed for the most inputs first
        rest = sylvan_project(r, game.ns_vars);
        std::vector<std::pair<double, int>> successors;
        std::vector<uint8_t> arr(game.statebits);
        for (MTBDD leaf = mtbdd_enum_all_first(rest, game.ns_vars, arr.data(), nullptr); leaf != mtbdd_false;
             leaf = mtbdd_enum_all_next(rest, game.ns_vars, arr.data(), nullptr)) {
            int c = 0;
            for (int i = 0; i < game.statebits; i++) c = 2 * c + (arr[i] ? 1 : 0);
            successors.emplace_back(0.0, index.at(c));
        }
        for (auto& s : successors) {
            lit = BDDTools::encode_state(code[s.second], game.ns_vars);
            can1 = sylvan_and_exists(r, lit, out);
            s.first = -sylvan_satcount(can1, game.uap_vars);
        }
        std::sort(successors.begin(), successors.end());
        for (const auto& s : successors) {
            const int t = s.second;
            lit = BDDTools::encode_state(code[t], game.ns_vars);
            can1 = sylvan_and_exists(r, lit, out);
            rest = sylvan_and(r, sylvan_not(lit));
            can0 = sylvan_exists(rest, out);
            const MTBDD g = choose(latches[t], k);
            link = sylvan_equiv(lit, g);
            r = sylvan_and(r, link);
        }
    };
    for (int k = 0; k < n; k++) {
        if (statesFirst) {
            chooseSuccessor(k);
            chooseOutputs(k);
        } else {
            chooseOutputs(k);
            chooseSuccessor(k);
        }
        // every input must still have a move (with one successor, as the latch functions select it)
        test = sylvan_exists(*rel[k], out);
        if (test != mtbdd_true) {
            std::cerr << "error: the one-hot encoding with don't cares lost the moves of a state" << std::endl;
            exit(1);
        }
    }

    // the circuit
    std::vector<unsigned int> latch(n), active(n);
    for (int k = 0; k < n; k++) {
        latch[k] = circuit->makeLatch();
        active[k] = code[k] == 0 ? aiger_not(latch[k]) : latch[k];
    }
    // a signal: the OR over its functions of "the function and one of the states that use it", where
    // with more than half of the states "none of the others" is smaller (exactly one state is active)
    auto signal = [&](const Functions& f) {
        std::deque<unsigned int> terms;
        for (size_t j = 0; j < f.function.size(); j++) {
            std::deque<unsigned int> sources;
            unsigned int guard;
            if (rules.guardComplement && 2 * f.states[j].size() > static_cast<size_t>(n)) {
                std::vector<char> uses(n, 0);
                for (int k : f.states[j]) uses[k] = 1;
                for (int k = 0; k < n; k++) {
                    if (!uses[k]) sources.push_back(active[k]);
                }
                guard = aiger_not(circuit->makeMultiOr(sources));
            } else {
                for (int k : f.states[j]) sources.push_back(active[k]);
                guard = circuit->makeMultiOr(sources);
            }
            const unsigned int inputs = isop ? bddToAigIsop(*f.function[j]) : bddToAigRecursive(*f.function[j]);
            terms.push_back(circuit->makeAnd(inputs, guard));
        }
        return circuit->makeMultiOr(terms);
    };
    for (int i = 0; i < game.cap_count; i++) circuit->makeOutput(signal(outputs[i]), cap_labels[i]);
    for (int t = 0; t < n; t++) {
        const unsigned int next = signal(latches[t]);
        circuit->setLatch(latch[t], code[t] == 0 ? aiger_not(next) : next, "");
    }
    return true;
}

namespace {

// The number of AND gates that bddToAigRecursive makes for the functions together: equal
// sub-BDDs give equal gates, as AIGCircuit hashes them. Adds the number of BDD nodes visited
// (a measure of the work, also of making the functions) to work.
class ShannonCost final {
public:
    static unsigned int of(const std::vector<MTBDD>& functions, long& work)
    {
        ShannonCost cost;
        for (MTBDD f : functions) cost.literal(f);
        work += static_cast<long>(cost.nodes.size()) + 1;
        return cost.gates;
    }

private:
    std::unordered_map<MTBDD, unsigned int> nodes;
    std::unordered_map<uint32_t, unsigned int> variables;
    std::unordered_map<uint64_t, unsigned int> ands;
    unsigned int next = 2;
    unsigned int gates = 0;

    unsigned int makeAnd(unsigned int a, unsigned int b)
    {
        if (b < a) std::swap(a, b);
        if (a == aiger_false) return aiger_false;
        if (a == aiger_true) return b;
        const uint64_t key = (static_cast<uint64_t>(b) << 32) | a;
        auto it = ands.find(key);
        if (it != ands.end()) return it->second;
        gates++;
        next += 2;
        return ands[key] = next - 2;
    }

    unsigned int literal(MTBDD bdd)
    {
        if (bdd == mtbdd_true) return aiger_true;
        if (bdd == mtbdd_false) return aiger_false;
        const bool complemented = (bdd & sylvan_complement) != 0;
        const MTBDD node = bdd & ~sylvan_complement;
        auto it = nodes.find(node);
        unsigned int res;
        if (it != nodes.end()) {
            res = it->second;
        } else {
            auto var = variables.find(mtbdd_getvar(node));
            if (var == variables.end()) {
                var = variables.emplace(mtbdd_getvar(node), next).first;
                next += 2;
            }
            const unsigned int v = var->second;
            const MTBDD low = mtbdd_getlow(node);
            const MTBDD high = mtbdd_gethigh(node);
            if (low == mtbdd_false) {
                res = makeAnd(v, literal(high));
            } else if (high == mtbdd_false) {
                res = makeAnd(aiger_not(v), literal(low));
            } else if (high == mtbdd_true) {
                res = aiger_not(makeAnd(aiger_not(v), aiger_not(literal(low))));
            } else {
                const unsigned int l = literal(low);
                const unsigned int h = literal(high);
                res = aiger_not(makeAnd(aiger_not(makeAnd(aiger_not(v), l)), aiger_not(makeAnd(v, h))));
            }
            nodes.emplace(node, res);
        }
        return complemented ? aiger_not(res) : res;
    }
};

// Choose a variable order for the functions, which the Shannon expansion follows. Sylvan has no
// dynamic reordering, so the functions are renamed to fresh variables (first, first + 1, ...) in
// the order to evaluate. Sifting moves each variable through all positions and keeps the best,
// until no move improves; then it restarts from a random order (with a fixed seed). The cost is
// the number of gates of the Shannon expansion of all functions together. It stops when the
// work (BDD nodes visited by the cost, see ShannonCost) exceeds budget, so the result does not
// depend on the speed of the machine, and adds the work to used. Returns false (and leaves the
// functions) if there are more than MAX_REORDER_VARS variables; otherwise replaces the
// functions and sets order[position] to the original variable.
constexpr int MAX_REORDER_VARS = 22;
constexpr int MAX_REORDER_RESTARTS = 3;

bool reorderFunctions(const std::vector<MTBDD*>& functions, const std::vector<uint32_t>& vars, uint32_t first,
                      long budget, long& used, std::vector<uint32_t>& order)
{
    const int n = static_cast<int>(vars.size());
    if (n == 0 || n > MAX_REORDER_VARS) return false;
    Pool pool;
    std::vector<MTBDD*> original, renamed;
    for (MTBDD* f : functions) {
        original.push_back(&pool.add(*f));
        renamed.push_back(&pool.add(mtbdd_false));
    }
    MTBDD& map = pool.add(mtbdd_false);
    MTBDD& var = pool.add(mtbdd_false);
    std::vector<MTBDD> values(functions.size());
    auto cost = [&](const std::vector<uint32_t>& ord) {
        map = mtbdd_map_empty();
        for (int pos = 0; pos < n; pos++) {
            var = sylvan_ithvar(first + static_cast<uint32_t>(pos));
            map = mtbdd_map_add(map, ord[pos], var);
        }
        for (size_t i = 0; i < functions.size(); i++) {
            *renamed[i] = sylvan_compose(*original[i], map);
            values[i] = *renamed[i];
        }
        return ShannonCost::of(values, used);
    };
    const long limit = used + budget;
    std::vector<uint32_t> current = vars;
    unsigned int currentCost = cost(current);
    order = current;
    unsigned int best = currentCost;
    uint32_t seed = 12345;
    auto random = [&](int bound) {
        seed = seed * 1103515245u + 12345u;
        return static_cast<int>((seed >> 16) % static_cast<uint32_t>(bound));
    };
    for (int restart = 0; restart < MAX_REORDER_RESTARTS; restart++) {
        if (restart > 0) {
            for (int i = n - 1; i > 0; i--) std::swap(current[i], current[random(i + 1)]);
            currentCost = cost(current);
        }
        bool improved = true;
        while (improved && used < limit) {
            improved = false;
            for (int v = 0; v < n && used < limit; v++) {
                std::vector<uint32_t> moved = current;
                moved.erase(std::find(moved.begin(), moved.end(), vars[v]));
                std::vector<uint32_t> bestMove = current;
                for (int pos = 0; pos < n; pos++) {
                    std::vector<uint32_t> candidate = moved;
                    candidate.insert(candidate.begin() + pos, vars[v]);
                    const unsigned int c = cost(candidate);
                    if (c < currentCost) {
                        currentCost = c;
                        bestMove = candidate;
                        improved = true;
                    }
                }
                current = bestMove;
            }
        }
        if (currentCost < best) {
            best = currentCost;
            order = current;
        }
        if (used >= limit) break;
    }
    cost(order);
    for (size_t i = 0; i < functions.size(); i++) *functions[i] = *renamed[i];
    return true;
}

}  // namespace

// After reorderFunctions (order[position] is the variable renamed to first +
// position), choose the functions again in that order: the choice among the
// simplifications of each bit (by node count) then sees the order of the
// circuit. State bits that are not live stay 0. Keeps the new functions if
// their Shannon expansion is smaller.
void AIGEncoder::rechoose(MTBDD relation, const std::vector<bool>& live, const std::vector<uint32_t>& order,
                          uint32_t first, MTBDD* cap_bdds, MTBDD* state_bdds)
{
    MTBDD map = mtbdd_map_empty(), var = mtbdd_false, renamed = mtbdd_false, inputs = mtbdd_set_empty();
    for (MTBDD* r : {&map, &var, &renamed, &inputs}) mtbdd_refs_pushptr(r);
    std::vector<uint32_t> inputPositions;
    const uint32_t firstUap = mtbdd_set_first(game.uap_vars);
    for (size_t pos = 0; pos < order.size(); pos++) {
        var = sylvan_ithvar(first + static_cast<uint32_t>(pos));
        map = mtbdd_map_add(map, order[pos], var);
        if (order[pos] >= firstUap && order[pos] < firstUap + static_cast<uint32_t>(game.uap_count)) {
            inputPositions.push_back(first + static_cast<uint32_t>(pos));
        }
    }
    for (int i = 0; i < game.statebits; i++) {
        if (!live[i]) map = mtbdd_map_add(map, mtbdd_set_first(game.s_vars) + static_cast<uint32_t>(i), mtbdd_false);
    }
    renamed = sylvan_compose(relation, map);
    for (int i = 0; i < game.statebits; i++) {
        if (live[i]) continue;
        var = sylvan_nithvar(mtbdd_set_first(game.ns_vars) + static_cast<uint32_t>(i));
        renamed = sylvan_and(renamed, var);
    }
    inputs = mtbdd_set_from_array(inputPositions.data(), inputPositions.size());
    std::vector<MTBDD> caps(game.cap_count), states(game.statebits);
    chooseWithDontCares(renamed, inputs, caps.data(), states.data());
    std::vector<MTBDD> before, after;
    for (int i = 0; i < game.cap_count; i++) {
        before.push_back(cap_bdds[i]);
        after.push_back(caps[i]);
    }
    for (int i = 0; i < game.statebits; i++) {
        if (!live[i]) continue;
        before.push_back(state_bdds[i]);
        after.push_back(states[i]);
    }
    if (ShannonCost::of(after, reorderUsed) < ShannonCost::of(before, reorderUsed)) {
        for (int i = 0; i < game.cap_count; i++) cap_bdds[i] = caps[i];
        for (int i = 0; i < game.statebits; i++) state_bdds[i] = states[i];
    }
    for (MTBDD& f : caps) mtbdd_unprotect(&f);
    for (MTBDD& f : states) mtbdd_unprotect(&f);
    mtbdd_refs_popptr(4);
}

void AIGEncoder::processBinary()
{
    MTBDD cap_bdds[game.cap_count];   // contains the solution: controllable ap bdds: state -> uap -> B
    MTBDD state_bdds[game.statebits]; // contains the solution: state bit bdds      : state -> uap -> B
    std::vector<bool> live(game.statebits, true);   // the state bits that get a latch
    std::vector<uint32_t> order;                     // with reordering: the variable at each position

    // we use the "strategies" field of the game class, which should only include state -> uap -> cap.

    MTBDD relation = mtbdd_false;   // with don't cares: the moves (state, input -> output, next state)
    mtbdd_refs_pushptr(&relation);
    if (dontcares) {
        relation = sylvan_and_exists(game.strategies, game.trans, game.np_vars);
        chooseWithDontCares(relation, game.uap_vars, cap_bdds, state_bdds);
        // A state bit whose next function is 0 is always 0 (latches start at 0): substitute 0 for
        // it in all functions and drop its latch, until no more bits become constant.
        MTBDD cube = mtbdd_false, set = mtbdd_false;
        mtbdd_refs_pushptr(&cube);
        mtbdd_refs_pushptr(&set);
        for (bool changed = rules.deadLatches; changed;) {
            changed = false;
            for (int i = 0; i < game.statebits; i++) {
                if (!live[i] || state_bdds[i] != mtbdd_false) continue;
                live[i] = false;
                changed = true;
                uint32_t var = mtbdd_set_first(game.s_vars) + static_cast<uint32_t>(i);
                cube = sylvan_nithvar(var);
                set = mtbdd_set_from_array(&var, 1);
                for (int j = 0; j < game.cap_count; j++) cap_bdds[j] = sylvan_and_exists(cap_bdds[j], cube, set);
                for (int j = 0; j < game.statebits; j++) {
                    if (live[j]) state_bdds[j] = sylvan_and_exists(state_bdds[j], cube, set);
                }
            }
        }
        mtbdd_refs_popptr(2);
        if (reorderWork > 0 && !isop) {
            std::vector<MTBDD*> functions;
            std::vector<uint32_t> vars;
            for (int i = 0; i < game.cap_count; i++) functions.push_back(&cap_bdds[i]);
            for (int i = 0; i < game.statebits; i++) {
                if (!live[i]) continue;
                functions.push_back(&state_bdds[i]);
                vars.push_back(mtbdd_set_first(game.s_vars) + static_cast<uint32_t>(i));
            }
            for (int i = 0; i < game.uap_count; i++) vars.push_back(mtbdd_set_first(game.uap_vars) + static_cast<uint32_t>(i));
            // fresh variables after all variables of the game (the next-state variables are last)
            const uint32_t first = mtbdd_set_first(game.ns_vars) + static_cast<uint32_t>(game.statebits);
            if (!reorderFunctions(functions, vars, first, reorderWork, reorderUsed, order)) {
                order.clear();
            } else if (rules.rechoose) {
                rechoose(relation, live, order, first, cap_bdds, state_bdds);
            }
        }
    } else {
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
    }

    // the latches of the state bits, and the literals of the variables (also renamed ones)
    {
        auto _vars = game.s_vars;
        for (int i = 0; i < game.statebits; i++) {
            auto bddvar = mtbdd_set_first(_vars);
            _vars = mtbdd_set_next(_vars);
            if (!live[i]) continue;
            auto lit = circuit->makeLatch();
            state_to_lit[i] = lit;
            bddvar_to_lit.emplace(bddvar, lit);
        }
        const uint32_t first = mtbdd_set_first(game.ns_vars) + static_cast<uint32_t>(game.statebits);
        for (size_t pos = 0; pos < order.size(); pos++) {
            bddvar_to_lit.emplace(first + static_cast<uint32_t>(pos), bddvar_to_lit.at(order[pos]));
        }
    }

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
            if (live[i]) {
                auto res = bddToAigCoverSop(state_zdds[i]);
                circuit->setLatch(state_to_lit.at(i), res, "");
            }
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
        std::vector<unsigned int> output_lits;
        for (int i=0; i<game.cap_count; i++) {
            auto res = bddToAigRecursive(cap_bdds[i]);
            circuit->makeOutput(res, cap_labels[i]);
            output_lits.push_back(res);
        }
        for (int i=0; i<game.statebits; i++) {
            if (!live[i]) continue;
            auto res = dontcares && rules.latchSplit ? latchFunction(state_bdds[i], cap_bdds, output_lits)
                                                     : bddToAigRecursive(state_bdds[i]);
            circuit->setLatch(state_to_lit.at(i), res, "");
        }
    }

    for (int i=0; i<game.cap_count; i++) mtbdd_unprotect(&cap_bdds[i]);
    for (int i=0; i<game.statebits; i++) mtbdd_unprotect(&state_bdds[i]);
    mtbdd_refs_popptr(1);
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
    } else if (onehot && dontcares) {
        if (!processOnehotDontCares(MAX_ONEHOT_STATES)) return nullptr;
    } else if (onehot) {
        processOnehot();
    } else {
        processBinary();
    }

    return std::move(circuit);
}
