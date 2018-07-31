#ifndef  JACOBIAN_FIELD
#define  JACOBIAN_FIELD

#include "TensorFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

class JacobianField: public TensorFiniteVolumeField
{
public:

    explicit JacobianField(const VectorFiniteVolumeField& u, const std::shared_ptr<CellGroup> &cells);

    void computeFaces();

    void compute(const CellGroup& cells);

    void compute();

private:

    const VectorFiniteVolumeField& u_;

};

#endif
