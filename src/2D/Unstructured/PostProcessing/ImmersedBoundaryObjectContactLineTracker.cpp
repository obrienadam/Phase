#include "FiniteVolumeGrid2D/BilinearInterpolator.h"

#include "ImmersedBoundaryObjectContactLineTracker.h"

ImmersedBoundaryObjectContactLineTracker::ImmersedBoundaryObjectContactLineTracker(int fileWriteFreq,
                                                                                   const std::weak_ptr<const ScalarFiniteVolumeField> &gamma,
                                                                                   const std::weak_ptr<const ImmersedBoundary> &ib)
    :
      Object(fileWriteFreq),
      gamma_(gamma),
      ib_(ib)
{
    path_ /= "ImmersedBoundaryObjectContactLineTracker";

    if (gamma_.lock()->grid()->comm().isMainProc())
        createOutputDirectory();

    for (const auto &ibObj: *ib_.lock())
    {
        if (gamma_.lock()->grid()->comm().isMainProc())
        {
            std::ofstream fout((path_ / (ibObj->name() + "_contact_lines.csv")).string());
            fout << "time,x,y,rx,ry,beta,nx,ny,theta\n";
            fout.close();
        }
    }
}

void ImmersedBoundaryObjectContactLineTracker::compute(Scalar time, bool force)
{
    typedef std::tuple<Point2D, Scalar> IbVal;

    if (do_update() || force)
    {
        for (const auto &ibObj: *ib_.lock())
        {
            std::vector<std::tuple<Point2D, Scalar>> ibVals;

            const auto &gamma = *gamma_.lock();
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
            ibVals = gamma_.lock()->grid()->comm().gatherv(gamma_.lock()->grid()->comm().mainProcNo(), ibVals);

            if (gamma_.lock()->grid()->comm().isMainProc())
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

                std::ofstream fout((path_ / (ibObj->name() + "_contact_lines.csv")).string(),
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
