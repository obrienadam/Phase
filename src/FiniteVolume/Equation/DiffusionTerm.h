#ifndef DIFFUSION_TERM_H
#define DIFFUSION_TERM_H

#include "Term.h"
#include "FiniteVolumeGrid2D.h"

class DiffusionTerm : public Term
{
public:
    DiffusionTerm(const ScalarFiniteVolumeField& var);
};

DiffusionTerm laplacian(const ScalarFiniteVolumeField& var);
DiffusionTerm laplacian(const ScalarFiniteVolumeField& coeffs, const ScalarFiniteVolumeField& var);

#endif
