#ifndef PHASE_POST_PROCESSING_H
#define PHASE_POST_PROCESSING_H

#include "Solvers/Solver.h"
#include "System/Input.h"
#include "System/PostProcessingInterface.h"

#include "CgnsViewer.h"

class PostProcessing : public PostProcessingInterface {
public:
  PostProcessing(const Input &input, const Solver &solver);

  // void initIbPostProcessingObjects(const Input &input, const Solver &solver);

  void compute(Scalar time, bool force = false) override;

protected:
  int _iter, _fileWriteFrequency;

  CgnsViewer _viewer;
};

#endif
