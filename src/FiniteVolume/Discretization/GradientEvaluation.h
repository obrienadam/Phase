#ifndef GRADIENT_EVALUATION_H
#define GRADIENT_EVALUATION_H

#include <functional>

#include "VectorFiniteVolumeField.h"
#include "TensorFiniteVolumeField.h"

namespace fv {

    enum GradientEvaluationMethod {
        GREEN_GAUSS_CELL_CENTERED, FACE_TO_CELL
    };

//- Scalar gradients
    void computeGradient(GradientEvaluationMethod method, const CellGroup &group, ScalarFiniteVolumeField &phi,
                         VectorFiniteVolumeField &gradPhi);

    void computeInverseWeightedGradient(const ScalarFiniteVolumeField &w, ScalarFiniteVolumeField &field,
                                        VectorFiniteVolumeField &gradField);

//- Field variant
    class Gradient: public VectorFiniteVolumeField
    {
    public:
        explicit Gradient(const ScalarFiniteVolumeField& phi);

        void computeFaces();

        void compute(const CellGroup& cells);

        void computeWeighted(const ScalarFiniteVolumeField& w, const CellGroup &group);

    private:
        const ScalarFiniteVolumeField& phi_;
    };
}

#endif
