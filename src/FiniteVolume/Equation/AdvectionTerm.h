#ifndef ADVECTION_TERM_H
#define ADVECTION_TERM_H

#include "Term.h"

class AdvectionTerm : public Term
{
public:
    AdvectionTerm(const ScalarFiniteVolumeField& var);
};

AdvectionTerm div(const ScalarFiniteVolumeField& var);

#endif
