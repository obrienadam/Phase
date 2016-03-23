#ifndef VIEWER_H
#define VIEWER_H

#include <vector>
#include <fstream>

#include "Solver.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

class Viewer
{
public:

    Viewer(const Solver& solver, const Input& input);
    ~Viewer();

    virtual void write(Scalar solutionTime) = 0;

    void addFieldToOutput(const ScalarFiniteVolumeField& field);
    void addFieldToOutput(const VectorFiniteVolumeField& field);

protected:

    const FiniteVolumeGrid2D& grid_;
    std::vector<const ScalarFiniteVolumeField*> scalarFields_;
    std::vector<const VectorFiniteVolumeField*> vectorFields_;

    std::string caseName_, outputFilename_;
    std::ofstream fout_;
};

#endif
