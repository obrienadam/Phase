#include "PostProcessing.h"
#include "CgnsViewer.h"
#include "CompactCgnsViewer.h"
#include "IbTracker.h"
#include "ImmersedBoundaryObjectContactLineTracker.h"
#include "ImmersedBoundaryObjectProbe.h"

PostProcessing::PostProcessing(const CommandLine &cl, const Input &input,
                               const Solver &solver) {
  iter_ = 0;
  fileWriteFrequency_ =
      input.postProcessingInput().get<int>("PostProcessing.fileWriteFrequency");

  std::string viewerType = input.postProcessingInput().get<std::string>(
      "PostProcessing.viewerType", "cgns");

  if (viewerType == "cgns")
    viewer_ = std::unique_ptr<Viewer>(new CgnsViewer(cl, input, solver));
  else if (viewerType == "compactCgns")
    viewer_ = std::unique_ptr<Viewer>(new CompactCgnsViewer(cl, input, solver));
  else
    throw Exception("PostProcessing", "PostProcessing",
                    "Unrecognized viewer type \"" + viewerType + "\".");
}

void PostProcessing::initIbPostProcessingObjects(const Input &input,
                                                 const Solver &solver) {
  auto objInputs =
      input.postProcessingInput().get_child_optional("PostProcessing.Objects");

  if (!objInputs || !solver.ib())
    return;

  for (const auto &objInput : objInputs.get()) {
    const std::string &name = objInput.first;
    const auto &inputTree = objInput.second;

    if (name == "IbTracker") {
      objs_.push_back(std::make_shared<IbTracker>(
          objInput.second.get<int>("fileWriteFrequency", fileWriteFrequency_),
          solver.ib()));
    } else if (name == "ImmersedBoundaryObjectProbe") {
      objs_.push_back(std::make_shared<ImmersedBoundaryObjectProbe>(
          objInput.second.get<int>("fileWriteFrequency", fileWriteFrequency_),
          solver.ib()->ibObj(inputTree.get<std::string>("name")),
          solver.scalarField(inputTree.get<std::string>("field")),
          inputTree.get<std::string>("position")));
    } else if (name == "ImmersedBoundaryObjectContactLineTracker") {
      objs_.push_back(
          std::make_shared<ImmersedBoundaryObjectContactLineTracker>(
              objInput.second.get<int>("fileWriteFrequency",
                                       fileWriteFrequency_),
              solver.scalarField(inputTree.get<std::string>("field", "gamma")),
              solver.ib()));
    }
  }
}

void PostProcessing::compute(Scalar time, bool force) {
  PostProcessingInterface::compute(time, force);

  if (iter_++ % fileWriteFrequency_ == 0 || force)
    viewer_->write(time);
}
