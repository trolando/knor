/**
 * Copyright Tom van Dijk
 */

#include <oink/oink.hpp>

extern "C" {
#include "simplehoa.h"
}

#ifndef KNOR_GAMECONSTRUCTOR_HPP
#define KNOR_GAMECONSTRUCTOR_HPP

class GameConstructor final {
public:
    /**
     * Construct a Oink parity game semi-symbolically.
     */
    static pg::Game* constructGame(HoaData* data, bool isMaxParity, bool controllerIsOdd);

    /**
     * Construct a Oink parity game fully explicitly.
     */
    static pg::Game* constructGameNaive(HoaData* data, bool isMaxParity, bool controllerIsOdd);

    static std::unique_ptr<SymGame> constructSymGame(HoaData* data, bool isMaxParity, bool controllerIsOdd);

    /**
     * This helper function ensures that the priority p is adjusted to
     * ensure we have a "max, even" parity game, as this is what Oink expects.
     * @param controllerIsOdd if the controller is the odd player
     * @param noPriorities how many priorities are in the game
     * @param maxPriority if the game is a max game
     */
    static int adjustPriority(int p, bool maxPriority, bool controllerIsOdd, int noPriorities);
};

#endif // KNOR_GAMECONSTRUCTOR_HPP
