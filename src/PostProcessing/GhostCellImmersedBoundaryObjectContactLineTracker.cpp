#include "GhostCellImmersedBoundaryObjectContactLineTracker.h"

GhostCellImmersedBoundaryObjectContactLineTracker::GhostCellImmersedBoundaryObjectContactLineTracker(
        const Solver &solver,
        std::shared_ptr<ImmersedBoundaryObject> gIbObj)
        :
        PostProcessing(solver)
{
    gIbObj_ = std::dynamic_pointer_cast<GhostCellImmersedBoundaryObject>(gIbObj);

    if (!gIbObj_)
    {
        solver.grid().comm().printf("WARNING: Immersed boundary object \"" + gIbObj->name() +
                                    "\" is not a ghost cell immersed boundary object.\n");
    }
    else
    {
        std::ofstream fout(gIbObj_->name() + "_contactLines.dat");
        fout << "time\tpos_x\tpos_y\tn_x\tn_y\ttheta\n";
        fout.close();
    }
}

void GhostCellImmersedBoundaryObjectContactLineTracker::compute(Scalar time)
{
    const ScalarFiniteVolumeField &gamma = solver_.scalarField("gamma");
    std::vector<std::pair<GhostCellStencil, Scalar>> gammaIb;

    for (const GhostCellStencil &st: gIbObj_->stencils())
        gammaIb.push_back(std::make_pair(
                st, st.bpValue(gamma)
        ));

    //- Sort ccw
    std::sort(gammaIb.begin(), gammaIb.end(), [this](const std::pair<GhostCellStencil, Scalar> &valA,
                                                     const std::pair<GhostCellStencil, Scalar> &valB) {
        Vector2D rA = valA.first.boundaryPoint() - gIbObj_->shape().centroid();
        Vector2D rB = valB.first.boundaryPoint() - gIbObj_->shape().centroid();
        Scalar thetaA = std::atan2(rA.y, rA.x);
        Scalar thetaB = std::atan2(rB.y, rB.x);
        return (thetaA < 0. ? thetaA + 2 * M_PI : thetaA) < (thetaB < 0. ? thetaB + 2 * M_PI : thetaB);
    });

    std::vector<std::pair<Point2D, Vector2D>> clData;
    for (int i = 0; i < gammaIb.size(); ++i)
    {
        const auto &ptA = gammaIb[i];
        const auto &ptB = gammaIb[(i + 1) % gammaIb.size()];
        bool isCandidate = ptA.second < 0.5 != ptB.second <= 0.5;

        if (isCandidate)
        {
            Scalar alpha = (0.5 - ptB.second) / (ptA.second - ptB.second);

            Point2D xc = gIbObj_->shape().nearestIntersect(
                    alpha * ptA.first.boundaryPoint()
                    + (1. - alpha) * ptB.first.boundaryPoint()
            );

            clData.push_back(std::make_pair(xc, Vector2D()));
        }
    }

    std::ofstream fout(gIbObj_->name() + "_contactLines.dat", std::ofstream::out | std::ofstream::app);

//    fout << "Zone T=\"Time=" << time << "s\"\n"
////         << "STRANDID=1, SOLUTIONTIME=" << time << "\n"
//         << "I = " << gammaIb.size() << "\n"
//         << "F=POINT\n";

    for (const auto &cl: clData)
    {
    }

    fout.close();
}
