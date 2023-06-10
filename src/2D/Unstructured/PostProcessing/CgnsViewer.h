#ifndef PHASE_CGNS_VIEWER_H
#define PHASE_CGNS_VIEWER_H

#include <vector>

#include "System/CgnsFile.h"

#include "Viewer.h"

class CgnsViewer : public Viewer {
public:
  CgnsViewer(const CommandLine &cl, const Input &input, const Solver &solver);

  virtual void write(Scalar time) override;

protected:
  std::string path_, gridfile_, casename_;
};

#endif
