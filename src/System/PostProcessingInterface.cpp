#include "PostProcessingInterface.h"

PostProcessingInterface::Object::Object()
{
    path_ = "./solution/PostProcessing";
}

bool PostProcessingInterface::Object::do_update()
{
    return iter_++ % fileWriteFreq_ == 0;
}

void PostProcessingInterface::Object::createOutputDirectory() const
{
    boost::filesystem::create_directories(path_);
}

void PostProcessingInterface::compute(Scalar time)
{
    for (auto &obj: objs_)
        obj->compute(time);
}