#include "FiniteVolume/Multiphase/CelesteImmersedBoundary.h"

#include "ImmersedBoundaryObjectContactLineTracker.h"

ImmersedBoundaryObjectContactLineTracker::
    ImmersedBoundaryObjectContactLineTracker(
        int fileWriteFreq,
        const std::weak_ptr<const ScalarFiniteVolumeField> &gamma,
        const std::weak_ptr<const ImmersedBoundary> &ib)
    : Object(fileWriteFreq), gamma_(gamma), ib_(ib) {
  path_ /= "ImmersedBoundaryObjectContactLineTracker";

  if (gamma_.lock()->grid()->comm().isMainProc())
    createOutputDirectory();

  for (const auto &ibObj : *ib_.lock()) {
    if (gamma_.lock()->grid()->comm().isMainProc()) {
      std::ofstream fout(
          (path_ / (ibObj->name() + "_contact_lines.csv")).string());
      fout << "time,x,y,rx,ry,beta,nx,ny,theta\n";
      fout.close();
    }
  }
}

void ImmersedBoundaryObjectContactLineTracker::compute(Scalar time,
                                                       bool force) {
  struct ContactLinePoint {
    Point2D pt;
    Scalar gamma;
    Vector2D n;
  };

  if (do_update() || force) {
    for (const auto &ibObj : *ib_.lock()) {
      std::vector<ContactLinePoint> clPts;
      clPts.reserve(ibObj->solidCells().size());

      const auto &gamma = *gamma_.lock();

      for (const Cell &cell : ibObj->solidCells()) {
        bool computeContactLine = false;

        for (const auto &nb : cell.neighbours())
          if (!ibObj->isInIb(nb.cell()))
            computeContactLine = true;

        if (!computeContactLine)
          continue;

        auto st = CelesteImmersedBoundary::ContactLineStencil(
            *ibObj, cell.centroid(), 111. * M_PI / 180., gamma);

        clPts.push_back(ContactLinePoint{st.cl()[1], st.gamma(), st.ncl()});
      }

      //- Gather all data to the main proc
      clPts = gamma_.lock()->grid()->comm().gatherv(
          gamma_.lock()->grid()->comm().mainProcNo(), clPts);

      if (gamma_.lock()->grid()->comm().isMainProc()) {
        //- Sort ccw
        std::sort(
            clPts.begin(), clPts.end(),
            [&ibObj](const ContactLinePoint &lhs, const ContactLinePoint &rhs) {
              return (lhs.pt - ibObj->shape().centroid()).angle() <
                     (rhs.pt - ibObj->shape().centroid()).angle();
            });

        std::vector<Point2D> clLocs;
        for (int i = 0; i < clPts.size(); ++i) {
          const auto &a = clPts[i];
          const auto &b = clPts[(i + 1) % clPts.size()];
          bool isCandidate = (a.gamma < 0.5) != (b.gamma <= 0.5);

          if (isCandidate) {
            Scalar alpha = (0.5 - b.gamma) / (a.gamma - b.gamma);

            Point2D xc = ibObj->shape().nearestIntersect(alpha * a.pt +
                                                         (1. - alpha) * b.pt);

            clLocs.push_back(xc);
          }
        }

        std::ofstream fout(
            (path_ / (ibObj->name() + "_contact_lines.csv")).string(),
            std::ofstream::out | std::ofstream::app);

        for (const auto &cl : clLocs) {
          Vector2D r = cl - ibObj->shape().centroid();

          fout << time << "," << cl.x << "," << cl.y << "," << r.x << "," << r.y
               << "," << std::atan2(r.y, r.x) << "," << 0. << "," << 0. << ","
               << 0. << "\n";
        }

        fout.close();
      }
    }
  }
}
