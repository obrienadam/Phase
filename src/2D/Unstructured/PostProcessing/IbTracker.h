#ifndef PHASE_IB_TRACKER_H
#define PHASE_IB_TRACKER_H

#include "PostProcessing.h"
#include "FiniteVolume/ImmersedBoundary/ImmersedBoundary.h"

class IbTracker : public PostProcessing::Object
{
public:

    IbTracker(int fileWriteFreq,
              const std::weak_ptr<const ImmersedBoundary> &ib,
              double lineThickness = 0.4,
              const std::string &fillColor = "CUST2");

    void compute(Scalar time, bool force = false) override;

private:

    std::weak_ptr<const ImmersedBoundary> ib_;

    int zoneNo_ = 1;

    double lineThickness_ = 0.4;

    std::string fillColor_ = "CUST2";
};


#endif
