#ifndef  JACOBIAN_FIELD
#define  JACOBIAN_FIELD

#include "TensorFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

class JacobianField: public TensorFiniteVolumeField
{
public:

    explicit JacobianField(const VectorFiniteVolumeField& u);

    void computeFaces();

    void compute(const CellGroup& cells);

private:

    const VectorFiniteVolumeField& u_;

};

#endif