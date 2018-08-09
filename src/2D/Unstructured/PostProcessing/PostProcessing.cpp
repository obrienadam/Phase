#include "PostProcessing.h"
#include "IbTracker.h"
#include "ImmersedBoundaryObjectProbe.h"
#include "ImmersedBoundaryObjectContactLineTracker.h"

PostProcessing::PostProcessing(const Input &input, const Solver &solver)
    :
      viewer_(input, solver)
{
    iter_ = 0;
    fileWriteFrequency_ = input.postProcessingInput().get<int>("PostProcessing.fileWriteFrequency");
}

void PostProcessing::initIbPostProcessingObjects(const Input &input, const Solver &solver)
{
    auto objInputs = input.postProcessingInput().get_child_optional("PostProcessing.Objects");

    if(!objInputs || !solver.ib())
        return;

    for (const auto &objInput: objInputs.get())
    {
        const std::string &name = objInput.first;
        const auto &inputTree = objInput.second;

        if (name == "IbTracker")
        {
            objs_.push_back(
                        std::make_shared<IbTracker>(objInput.second.get<int>("fileWriteFrequency", fileWriteFrequency_), solver.ib())
                        );
        }
        else if (name == "ImmersedBoundaryObjectProbe")
        {
            objs_.push_back(
                        std::make_shared<ImmersedBoundaryObjectProbe>(
                            objInput.second.get<int>("fileWriteFrequency", fileWriteFrequency_),
                            solver.ib()->ibObj(inputTree.get<std::string>("name")),
                            solver.scalarField(inputTree.get<std::string>("field")),
                            inputTree.get<std::string>("position")
                            ));
        }
        else if (name == "ImmersedBoundaryObjectContactLineTracker")
        {
            objs_.push_back(
                        std::make_shared<ImmersedBoundaryObjectContactLineTracker>(
                            objInput.second.get<int>("fileWriteFrequency", fileWriteFrequency_),
                            solver.scalarField(inputTree.get<std::string>("field")),
                            solver.ib()
                            ));
        }
    }
}

void PostProcessing::compute(Scalar time, bool force)
{
    PostProcessingInterface::compute(time, force);

    if (++iter_ % fileWriteFrequency_ == 0 || force)
        viewer_.write(time);
}
