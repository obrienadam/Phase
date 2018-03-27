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

    auto objInputs = input.postProcessingInput().get_child_optional("PostProcessing.Objects");

    if(!objInputs)
        return;

    for (const auto &objInput: objInputs.get())
    {
        const std::string &name = objInput.first;
        const auto &inputTree = objInput.second;

        if (name == "IbTracker")
        {
            objs_.push_back(
                    std::make_shared<IbTracker>(solver)
            );
        }
        else if (name == "ImmersedBoundaryObjectProbe")
        {
            objs_.push_back(
                    std::make_shared<ImmersedBoundaryObjectProbe>(
                            solver,
                            inputTree.get<std::string>("name"),
                            inputTree.get<std::string>("field"),
                            inputTree.get<std::string>("position")
                    ));
        }
        else if (name == "ImmersedBoundaryObjectContactLineTracker")
        {
            objs_.push_back(
                    std::make_shared<ImmersedBoundaryObjectContactLineTracker>(solver)
            );
        }
    }
}

void PostProcessing::compute(Scalar time)
{
    PostProcessingInterface::compute(time);

    if (++iter_ % fileWriteFrequency_ == 0)
        viewer_.write(time);
}