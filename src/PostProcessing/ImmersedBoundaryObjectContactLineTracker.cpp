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
                std::ofstream fout((outputDir_ / (ibObj->name() + "_contact_lines.csv")).string());
                fout << "time,x,y,rx,ry,beta,nx,ny,theta\n";
                fout.close();
            }
        }
    }
}

void ImmersedBoundaryObjectContactLineTracker::compute(Scalar time)
{
    typedef std::tuple<Point2D, Scalar> IbVal;

    if (iterNo_++ % fileWriteFrequency_ == 0)
    {
        for (const auto &ibObj: solver_.ib())
        {
            std::vector<std::tuple<Point2D, Scalar>> ibVals;

            const auto &gamma = solver_.scalarField("gamma");
            auto bi = BilinearInterpolator(gamma.grid());
            for (const Cell &cell: ibObj->ibCells())
            {
                bi.setPoint(ibObj->nearestIntersect(cell.centroid()));

                if (bi.isValid())
                {
                    ibVals.push_back(std::make_tuple(
                            bi.point(),
                            bi(gamma)
                    ));
                }
            }

            //- Gather all data to the main proc
            ibVals = solver_.grid()->comm().gatherv(solver_.grid()->comm().mainProcNo(), ibVals);

            if (solver_.grid()->comm().isMainProc())
            {
                //- Sort ccw
                std::sort(ibVals.begin(), ibVals.end(), [&ibObj](const IbVal &lhs,
                                                                 const IbVal &rhs)
                {
                    return (std::get<0>(lhs) - ibObj->shape().centroid()).angle() <
                           (std::get<0>(rhs) - ibObj->shape().centroid()).angle();
                });

                std::vector<Vector2D> clPoints;
                for (int i = 0; i < ibVals.size(); ++i)
                {
                    const auto &a = ibVals[i];
                    const auto &b = ibVals[(i + 1) % ibVals.size()];
                    bool isCandidate = (std::get<1>(a) < 0.5) != (std::get<1>(b) < 0.5);

                    if (isCandidate)
                    {
                        Scalar alpha = (0.5 - std::get<1>(b)) / (std::get<1>(a) - std::get<1>(b));

                        Point2D xc = ibObj->shape().nearestIntersect(
                                alpha * std::get<0>(a) + (1. - alpha) * std::get<0>(b)
                        );

                        clPoints.push_back(xc);
                    }
                }

                std::ofstream fout((outputDir_ / (ibObj->name() + "_contact_lines.csv")).string(),
                                   std::ofstream::out | std::ofstream::app);

                for (const auto &cl: clPoints)
                {
                    Vector2D r = cl - ibObj->shape().centroid();

                    fout << time << ","
                         << cl.x << ","
                         << cl.y << ","
                         << r.x << ","
                         << r.y << ","
                         << std::atan2(r.y, r.x) << ","
                         << 0. << ","
                         << 0. << ","
                         << 0. << "\n";
                }

                fout.close();
            }
        }
    }
}
