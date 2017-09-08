#include "PostProcessingObject.h"

PostProcessingObject::PostProcessingObject(const Solver &solver)
        :
        solver_(solver)
{
    outputDir_ = "./solution/PostProcessing";
}

void PostProcessingObject::createOutputDirectory() const
{
    boost::filesystem::create_directories(outputDir_);
}