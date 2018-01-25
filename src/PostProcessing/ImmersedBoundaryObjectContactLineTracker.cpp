#include "ImmersedBoundaryObjectContactLineTracker.h"
#include "BilinearInterpolator.h"

ImmersedBoundaryObjectContactLineTracker::ImmersedBoundaryObjectContactLineTracker(const Solver &solver)
        :
        PostProcessingObject(solver)
{
    outputDir_ /= "ImmersedBoundaryObjectContactLineTracker";

    if (solver.grid()->comm().isMainProc())
        createOutputDirectory();

    for (const auto &ibObj: solver_.ib())
    {
        if (solver_.scalarFieldPtr("gamma"))
        {
            if (solver.grid()->comm().isMainProc())
            {
                std::ofstream fout((outputDir_ / (ibObj->name() + "_contactLines.dat")).string());
                fout << "time\tpos_x\tpos_y\tn_x\tn_y\ttheta\n";
                fout.close();
            }
        }
    }
}

void ImmersedBoundaryObjectContactLineTracker::compute(Scalar time)
{
    if (iterNo_++ % fileWriteFrequency_ == 0)
    {
        for (const auto &ibObj: solver_.ib())
        {
            std::vector<Vector2D> pts;
            std::vector<Scalar> gammaVals;

            const auto &gamma = solver_.scalarField("gamma");
            auto bi = BilinearInterpolator(gamma.grid());
            for (const Cell &cell: ibObj->ibCells())
            {
                bi.setPoint(ibObj->nearestIntersect(cell.centroid()));

                if (bi.isValid())
                {
                    pts.push_back(bi.point());
                    gammaVals.push_back(bi(gamma));
                }
            }

            //- Gather all data to the main proc
            pts = solver_.grid()->comm().gatherv(solver_.grid()->comm().mainProcNo(), pts);
            gammaVals = solver_.grid()->comm().gatherv(solver_.grid()->comm().mainProcNo(), gammaVals);

            if (solver_.grid()->comm().isMainProc())
            {
                std::vector<std::pair<Point2D, Scalar>> gammaBps(pts.size());
                std::transform(pts.begin(), pts.end(), gammaVals.begin(), gammaBps.begin(),
                               [](const Point2D &pt, Scalar val) { return std::make_pair(pt, val); });

                //- Sort ccw
                std::sort(gammaBps.begin(), gammaBps.end(), [&ibObj](const std::pair<Point2D, Scalar> &lhs,
                                                                     const std::pair<Point2D, Scalar> &rhs) {
                    return (lhs.first - ibObj->shape().centroid()).angle() <
                           (rhs.first - ibObj->shape().centroid()).angle();
                });

                std::vector<Vector2D> clPoints;
                for (int i = 0; i < gammaBps.size(); ++i)
                {
                    const auto &ptA = gammaBps[i];
                    const auto &ptB = gammaBps[(i + 1) % gammaBps.size()];
                    bool isCandidate = (ptA.second < 0.5) != (ptB.second < 0.5);

                    if (isCandidate)
                    {
                        Scalar alpha = (0.5 - ptB.second) / (ptA.second - ptB.second);

                        Point2D xc = ibObj->shape().nearestIntersect(
                                alpha * ptA.first + (1. - alpha) * ptB.first
                        );

                        clPoints.push_back(xc);
                    }
                }

                std::ofstream fout((outputDir_ / (ibObj->name() + "_contactLines.dat")).string(),
                                   std::ofstream::out | std::ofstream::app);

                for (const auto &cl: clPoints)
                {
                    Vector2D r = cl - ibObj->shape().centroid();
                    fout << time << "\t" << cl.x << "\t" << cl.y << "\t" << 0 << "\t" << 0 << "\t" << std::atan2(r.y, r.x) << "\n";
                }

                fout.close();
            }
        }
    }
}
