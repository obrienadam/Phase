#ifndef CONTACT_LINE_TRACKER_H
#define CONTACT_LINE_TRACKER_H

#include "Patch.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

class ContactLineTracker
{

    void compute(Scalar time);

    void write() const;

private:

    const Patch& finpatch_;

    std::vector<Scalar> time_;
    std::vector<Point2D> location_;
    std::vector<Scalar> angle_;
};

#endif