#include "GhostCellImmersedBoundaryObjectContactLineTracker.h"

GhostCellImmersedBoundaryObjectContactLineTracker::GhostCellImmersedBoundaryObjectContactLineTracker(
        const Solver &solver)
        :
        PostProcessingObject(solver)
{
    outputDir_ /= "GhostCellImmersedBoundaryObjectContactLineTracker";

    if (solver.grid().comm().isMainProc())
        createOutputDirectory();

    for (auto ibObj: solver_.ib().ibObjPtrs())
    {
        auto gcIbObj = std::dynamic_pointer_cast<GhostCellImmersedBoundaryObject>(ibObj);

        if (gcIbObj && solver_.scalarFieldPtr("gamma"))
        {
            if (solver.grid().comm().isMainProc())
            {
                std::ofstream fout((outputDir_ / (gcIbObj->name() + "_contactLines.dat")).string());
                fout << "time\tpos_x\tpos_y\tn_x\tn_y\ttheta\n";
                fout.close();
            }

            gcIbObjs_.push_back(gcIbObj);
        }
    }
}

void GhostCellImmersedBoundaryObjectContactLineTracker::compute(Scalar time)
{
    const ScalarFiniteVolumeField &gamma = solver_.scalarField("gamma");

    gcIbObjs_.erase(
            std::remove_if(gcIbObjs_.begin(), gcIbObjs_.end(),
                           [](std::weak_ptr<GhostCellImmersedBoundaryObject> &ibObj) {
                               return !ibObj.lock();
                           }), gcIbObjs_.end());

    if (iterNo_++ % fileWriteFrequency_ == 0)
    {
        for (auto ptr: gcIbObjs_)
        {
            auto gcIbObj = ptr.lock();

            std::vector<std::pair<GhostCellStencil, Scalar>> gammaIb;

            for (const GhostCellStencil &st: gcIbObj->stencils())
                gammaIb.push_back(std::make_pair(
                        st, st.bpValue(gamma)
                ));

            //- Sort ccw
            std::sort(gammaIb.begin(), gammaIb.end(), [gcIbObj](const std::pair<GhostCellStencil, Scalar> &valA,
                                                                const std::pair<GhostCellStencil, Scalar> &valB) {
                Vector2D rA = valA.first.boundaryPoint() - gcIbObj->shape().centroid();
                Vector2D rB = valB.first.boundaryPoint() - gcIbObj->shape().centroid();
                Scalar thetaA = std::atan2(rA.y, rA.x);
                Scalar thetaB = std::atan2(rB.y, rB.x);
                return (thetaA < 0. ? thetaA + 2 * M_PI : thetaA) < (thetaB < 0. ? thetaB + 2 * M_PI : thetaB);
            });

            std::vector<Vector2D> clPoints;
            for (int i = 0; i < gammaIb.size(); ++i)
            {
                if(i == gammaIb.size() - 1)
                    continue;

                const auto &ptA = gammaIb[i];
                const auto &ptB = gammaIb[(i + 1) % gammaIb.size()];
                bool isCandidate = ptA.second < 0.5 != ptB.second <= 0.5;

                if (isCandidate)
                {
                    Scalar alpha = (0.5 - ptB.second) / (ptA.second - ptB.second);

                    Point2D xc = gcIbObj->shape().nearestIntersect(
                            alpha * ptA.first.boundaryPoint()
                            + (1. - alpha) * ptB.first.boundaryPoint()
                    );

                    clPoints.push_back(xc);
                }
            }

            if (solver_.grid().comm().isMainProc())
            {

                std::ofstream fout((outputDir_ / (gcIbObj->name() + "_contactLines.dat")).string(),
                                   std::ofstream::out | std::ofstream::app);

                for (const auto &cl: clPoints)
                    fout << time << "\t" << cl.x << "\t" << cl.y << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\n";

                //- Probe other processes
                for(int proc = 0; proc < solver_.grid().comm().nProcs(); ++proc)
                {
                    if(proc == solver_.grid().comm().rank())
                        continue;

                    clPoints.resize(solver_.grid().comm().probeSize<Vector2D>(proc, proc));
                    solver_.grid().comm().recv(proc, clPoints, proc);

                    for (const auto &cl: clPoints)
                        fout << time << "\t" << cl.x << "\t" << cl.y << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\n";
                }

                fout.close();
            }
            else
            {
                solver_.grid().comm().ssend(solver_.grid().comm().mainProcNo(), clPoints, solver_.grid().comm().rank());
            }
        }
    }
}
