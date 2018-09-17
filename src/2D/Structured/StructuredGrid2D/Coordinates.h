#ifndef PHASE_COORDINATES_H
#define PHASE_COORDINATES_H

#include <array>

class Coordinates
{
public:

    enum Index{I, J};

    enum Direction{I_POS, I_NEG, J_NEG, J_POS};

    static std::array<Direction, 4> DIRECTIONS;
};

#endif
