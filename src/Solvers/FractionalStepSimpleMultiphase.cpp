#include "FractionalStepSimpleMultiphase.h"

FractionalStepSimpleMultiphase::FractionalStepSimpleMultiphase(const Input &input,
                                                               const Communicator &comm,
                                                               FiniteVolumeGrid2D &grid)
        :
        FractionalStepSimple(input,
                             comm,
                             grid),
        rho(addScalarField("rho")),
        mu(addScalarField("mu")),
        gamma(addScalarField(input, "gamma")),
        rhoU(addVectorField("rhoU")),
        gammaEqn_(input, comm, gamma, "gammaEqn")
{

}
