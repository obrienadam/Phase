#include "GhostCellImmersedBoundaryObjectContactLineTracker.h"

GhostCellImmersedBoundaryObjectContactLineTracker::GhostCellImmersedBoundaryObjectContactLineTracker(
        const Solver &solver)
        :
        PostProcessingObject(solver)
{
    outputDir_ /= "GhostCellImmersedBoundaryObjectContactLineTracker";

    if (solver.grid().comm().isMainProc())
        createOutputDirectory();

    for (const auto &ibObj: solver_.ib())
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
                           [](std::weak_ptr <GhostCellImmersedBoundaryObject> &ibObj) {
                               return !ibObj.lock();
                           }), gcIbObjs_.end());

    if (iterNo_++ % fileWriteFrequency_ == 0)
    {
        for (auto ptr: gcIbObjs_)
        {
            auto gcIbObj = ptr.lock();

            std::vector <Vector2D> pts;
            std::vector <Scalar> gammaVals;

            for (const GhostCellStencil &st: gcIbObj->stencils())
            {
                pts.push_back(st.boundaryPoint());
                gammaVals.push_back(st.bpValue(gamma));
            }

            //- Gather all data to the main proc
            pts = solver_.grid().comm().gatherv(solver_.grid().comm().mainProcNo(), pts);
            gammaVals = solver_.grid().comm().gatherv(solver_.grid().comm().mainProcNo(), gammaVals);

            if (solver_.grid().comm().isMainProc())
            {
                std::vector <std::pair<Point2D, Scalar>> gammaBps;
                std::transform(pts.begin(), pts.end(), gammaVals.begin(), std::back_inserter(gammaBps),
                               [](const Point2D &pt, Scalar val) { return std::make_pair(pt, val); });

                //- Sort ccw
                std::sort(gammaBps.begin(), gammaBps.end(), [gcIbObj](const std::pair <Point2D, Scalar> &valA,
                                                              const std::pair <Point2D, Scalar> &valB) {
                    Vector2D rA = valA.first - gcIbObj->shape().centroid();
                    Vector2D rB = valB.first - gcIbObj->shape().centroid();
                    return rA.angle() < rB.angle();
                });

                std::vector <Vector2D> clPoints;
                for (int i = 0; i < gammaBps.size(); ++i)
                {
                    const auto &ptA = gammaBps[i];
                    const auto &ptB = gammaBps[(i + 1) % gammaBps.size()];
                    bool isCandidate = (ptA.second < 0.5) != (ptB.second < 0.5);

                    if (isCandidate)
                    {
                        Scalar alpha = (0.5 - ptB.second) / (ptA.second - ptB.second);

                        Point2D xc = gcIbObj->shape().nearestIntersect(
                                alpha * ptA.first + (1. - alpha) * ptB.first
                        );

                        clPoints.push_back(xc);
                    }
                }

                std::ofstream fout((outputDir_ / (gcIbObj->name() + "_contactLines.dat")).string(),
                                   std::ofstream::out | std::ofstream::app);

                for (const auto &cl: clPoints)
                    fout << time << "\t" << cl.x << "\t" << cl.y << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\n";

                fout.close();
            }
        }
    }
}
