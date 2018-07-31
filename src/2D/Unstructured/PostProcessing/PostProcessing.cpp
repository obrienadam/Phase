#include "PostProcessing.h"
#include "IbTracker.h"
#include "ImmersedBoundaryObjectProbe.h"
#include "ImmersedBoundaryObjectContactLineTracker.h"

PostProcessing::PostProcessing(const Input &input, const Solver &solver, const std::weak_ptr<const ImmersedBoundary> &ib)
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

        if (name == "IbTracker" && ib.lock())
        {
            objs_.push_back(
                        std::make_shared<IbTracker>(objInput.second.get<int>("fileWriteFrequency", fileWriteFrequency_), ib)
                        );
        }
        else if (name == "ImmersedBoundaryObjectProbe" && ib.lock())
        {
            objs_.push_back(
                        std::make_shared<ImmersedBoundaryObjectProbe>(
                            objInput.second.get<int>("fileWriteFrequency", fileWriteFrequency_),
                            ib.lock()->ibObj(inputTree.get<std::string>("name")),
                            solver.scalarField(inputTree.get<std::string>("field")),
                            inputTree.get<std::string>("position")
                            ));
        }
        else if (name == "ImmersedBoundaryObjectContactLineTracker" && ib.lock())
        {
            objs_.push_back(
                        std::make_shared<ImmersedBoundaryObjectContactLineTracker>(
                            objInput.second.get<int>("fileWriteFrequency", fileWriteFrequency_),
                            solver.scalarField(inputTree.get<std::string>("field")),
                            ib
                            ));
        }
    }
}

void PostProcessing::compute(Scalar time)
{
    PostProcessingInterface::compute(time);

    if (++iter_ % fileWriteFrequency_ == 0)
        viewer_.write(time);
}
