#include "PostProcessing.h"
#include "VolumeIntegrator.h"
#include "IbTracker.h"
#include "GhostCellImmersedBoundaryObjectForceIntegrator.h"
#include "GhostCellImmersedBoundaryObjectContactLineTracker.h"

PostProcessing::PostProcessing(const Input &input, const Solver &solver)
{
    auto postProcessingInput = input.postProcessingInput().get_child_optional("PostProcessing");
    int defaultFileWriteFrequency = input.caseInput().get<int>("System.fileWriteFrequency", 1);

    if (postProcessingInput)
    {
        for (const auto &postProcessingObjectInput: postProcessingInput.get())
        {
            std::shared_ptr<PostProcessingObject> postProcessingObj;

            if (postProcessingObjectInput.first == "VolumeIntegrator")
            {
                postProcessingObj = std::make_shared<VolumeIntegrator>(
                        solver,
                        postProcessingObjectInput.second.get<std::string>("field")
                );
            }
            else if (postProcessingObjectInput.first == "IbTracker")
            {
                postProcessingObj = std::make_shared<IbTracker>(
                        solver,
                        postProcessingObjectInput.second.get<double>("lineThickness", 0.4),
                        postProcessingObjectInput.second.get<std::string>("fillColor", "CUST2")
                );
            }
            else if(postProcessingObjectInput.first == "GhostCellImmersedBoundaryObjectForceIntegrator")
            {
                postProcessingObj = std::make_shared<GhostCellImmersedBoundaryObjectForceIntegrator>(solver);
            }
            else if (postProcessingObjectInput.first == "GhostCellImmersedBoundaryObjectContactLineTracker")
            {
                postProcessingObj = std::make_shared<GhostCellImmersedBoundaryObjectContactLineTracker>(solver);
            }

            if (postProcessingObj)
            {
                postProcessingObjs_.push_back(postProcessingObj);
                postProcessingObj->setFileWriteFrequency(
                        postProcessingObjectInput.second.get<int>("fileWriteFrequency", defaultFileWriteFrequency)
                );
            }
        }
    }
}

void PostProcessing::compute(Scalar time) const
{
    for (auto obj: postProcessingObjs_)
        obj->compute(time);
}