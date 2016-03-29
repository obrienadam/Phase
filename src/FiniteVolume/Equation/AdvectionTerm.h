#ifndef ADVECTION_TERM_H
#define ADVECTION_TERM_H

#include "Term.h"

class AdvectionTerm : public Term
{
public:
    AdvectionTerm(const VectorFiniteVolumeField& u, const ScalarFiniteVolumeField& var);
    AdvectionTerm(const VectorFiniteVolumeField &u, const VectorFiniteVolumeField &var);
};

AdvectionTerm div(const VectorFiniteVolumeField& u, const ScalarFiniteVolumeField& var);
AdvectionTerm div(const VectorFiniteVolumeField& u, const VectorFiniteVolumeField& var);

#endif
