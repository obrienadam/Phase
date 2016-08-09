#ifndef VIEWER_H
#define VIEWER_H

#include <vector>

#include "Solver.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

class Viewer
{
public:

    Viewer(const Solver& solver, const Input& input);

    virtual void write(Scalar solutionTime);

protected:

    const Solver& solver_;

    int fileId_, baseId_, zoneId_;
    std::vector<Scalar> timeValues_;
    std::vector<std::string> flowSolutionPointers_;

    std::vector< Ref<const ScalarFiniteVolumeField> > scalarFields_;
    std::vector< Ref<const VectorFiniteVolumeField> > vectorFields_;

    std::string caseName_, outputFilename_;
};

#endif
