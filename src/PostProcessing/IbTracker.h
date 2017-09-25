#ifndef IB_TRACKER_H
#define IB_TRACKER_H

#include "PostProcessingObject.h"
#include "ImmersedBoundary.h"

class IbTracker : public PostProcessingObject
{
public:
    IbTracker(const Solver &solver,
              double lineThickness = 0.4,
              const std::string &fillColor = "CUST2");

    void compute(Scalar time);

private:
    std::vector<std::weak_ptr<ImmersedBoundaryObject>> ibObjs_;
    int zoneNo_ = 1;
    double lineThickness_ = 0.4;
    std::string fillColor_ = "CUST2";
};


#endif
