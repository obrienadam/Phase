#ifndef VIEWER_H
#define VIEWER_H

#include <vector>

#include "Solver.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

class Viewer
{
public:

    Viewer(const Input& input, const Solver& solver);
    Viewer(const Input& input, const Communicator& comm, const Solver& solver);

    virtual void write(Scalar solutionTime);
    virtual void write(Scalar solutionTime, const Communicator& comm);

    virtual void write(const std::vector<VolumeIntegrator>& volumeIntegrators);

protected:

    void init(const Input& input, const Solver& solver, const std::string& filename);
    void init(const Input& input, const Communicator& comm, const Solver& solver, const std::string& filename);

    void initFields(const Input& input);
    int  createBase(const std::string& name = "Case");
    int  createZone(const FiniteVolumeGrid2D& grid, const std::string& name = "Cells");

    void writeCoords(const FiniteVolumeGrid2D& grid);
    int  writeConnectivity(const FiniteVolumeGrid2D& grid);
    void writeBoundaryConnectivity(const FiniteVolumeGrid2D& grid);
    void writeImmersedBoundaries(const Solver& solver);

    void updateBaseIterativeData(Scalar solutionTime);
    int  updateFlowSolutionPointers();

    const Solver& solver_;

    std::string filename_;
    int fileId_, baseId_, zoneId_;
    std::vector<Scalar> timeValues_;
    std::vector<std::string> flowSolutionPointers_;

    std::vector< Ref<const ScalarFiniteVolumeField> > scalarFields_;
    std::vector< Ref<const VectorFiniteVolumeField> > vectorFields_;
};

#endif
