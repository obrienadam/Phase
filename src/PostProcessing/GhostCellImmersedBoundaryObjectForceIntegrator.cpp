#include "GhostCellImmersedBoundaryObjectForceIntegrator.h"

GhostCellImmersedBoundaryObjectForceIntegrator::GhostCellImmersedBoundaryObjectForceIntegrator(const Solver &solver)
        :
        PostProcessingObject(solver)
{
    outputDir_ /= "GhostCellImmersedBoundaryObjectForceIntegrator";

    if (solver.grid()->comm().isMainProc())
        createOutputDirectory();

    for (const auto &ibObj: solver_.ib())
    {
        auto gcIbObj = std::dynamic_pointer_cast<GhostCellImmersedBoundaryObject>(ibObj);

        if (gcIbObj)
        {
            if (solver.grid()->comm().isMainProc())
            {
                std::ofstream fout((outputDir_ / (gcIbObj->name() + "_contactLines.dat")).string());
                fout << "time\tf_x\tf_y\n";
                fout.close();
            }

            gcIbObjs_.push_back(gcIbObj);
        }
    }
}
void GhostCellImmersedBoundaryObjectForceIntegrator::compute(Scalar time)
{
    gcIbObjs_.erase(
            std::remove_if(gcIbObjs_.begin(), gcIbObjs_.end(), [](std::weak_ptr<GhostCellImmersedBoundaryObject> gcIbObj){
                return !gcIbObj.lock();
            }),
            gcIbObjs_.end()
    );

    for(auto ptr: gcIbObjs_)
    {
        auto gcIbObj = ptr.lock();


    }
}
