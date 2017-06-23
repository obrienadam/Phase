#ifndef VIEWER_H
#define VIEWER_H

#include "Input.h"
#include "Communicator.h"
#include "Solver.h"

class Viewer
{
public:

    Viewer(const Input& input, const Solver& solver);

    virtual void write(Scalar solutionTime) = 0;

protected:

    const Solver& solver_;

    std::string filename_;

    std::vector< Ref<const FiniteVolumeField<int> > > integerFields_;
    std::vector< Ref<const ScalarFiniteVolumeField> > scalarFields_;
    std::vector< Ref<const VectorFiniteVolumeField> > vectorFields_;

};

#include "CgnsViewer.h"
#include "IbViewer.h"

#endif
