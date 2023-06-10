#ifndef PHASE_VIEWER_H
#define PHASE_VIEWER_H

#include "Solvers/Solver.h"
#include "System/CommandLine.h"
#include "System/Communicator.h"
#include "System/Input.h"

class Viewer {
public:
  Viewer(const CommandLine &cl, const Input &input, const Solver &solver);

  virtual void write(Scalar solutionTime) = 0;

protected:
  const Solver &solver_;

  std::string filename_;

  std::unordered_set<std::string> integerFields_, scalarFields_, vectorFields_;

  bool isRestart_;
};

#include "CgnsViewer.h"

#endif
