#ifndef PISO_H
#define PISO_H

#include "Simple.h"

class Piso : public Simple
{
public:

    Piso(const FiniteVolumeGrid2D& grid, const Input& input);

    Scalar solve(Scalar timeStep);

private:

    size_t nPCorrections_;
};

#endif
