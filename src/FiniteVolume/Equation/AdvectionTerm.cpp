#include "AdvectionTerm.h"

AdvectionTerm::AdvectionTerm(const ScalarFiniteVolumeField &var)
    :
      Term(var)
{

}

//- External functions

AdvectionTerm div(const ScalarFiniteVolumeField &var)
{
    return AdvectionTerm(var);
}
