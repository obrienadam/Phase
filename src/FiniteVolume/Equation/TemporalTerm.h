#ifndef TEMPORAL_TERM_H
#define TEMPORAL_TERM_H

#include "Term.h"

class TemporalTerm : public Term
{
public:

    TemporalTerm(const ScalarFiniteVolumeField& var, Scalar timeStep);
    TemporalTerm(const VectorFiniteVolumeField& var, Scalar timeStep);

private:
};

TemporalTerm ddt(const ScalarFiniteVolumeField& var, Scalar timeStep);
TemporalTerm ddt(const VectorFiniteVolumeField& var, Scalar timeStep);

#endif
