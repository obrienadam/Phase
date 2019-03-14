#ifndef PHASE_COMPACT_CGNS_VIEWER_H
#define PHASE_COMPACT_CGNS_VIEWER_H

#include "Viewer.h"

class CompactCgnsViewer: public Viewer
{
public:

    CompactCgnsViewer(const CommandLine &cl, const Input& input, const Solver& solver);

    virtual void write(Scalar time) override;

protected:

    int bid_, zid_;

    std::size_t solnNo_;

    std::string filename_;

};

#endif
