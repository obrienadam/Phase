#ifndef POST_PROCESSING_OBJECT_H
#define POST_PROCESSING_OBJECT_H

#include <boost/filesystem.hpp>

#include "Solver.h"

class PostProcessingObject
{
public:

    PostProcessingObject(const Solver &solver);

    virtual void compute(Scalar time) = 0;

    void setFileWriteFrequency(int fileWriteFrequency)
    { fileWriteFrequency_ = fileWriteFrequency; }

protected:

    void createOutputDirectory() const;

    const Solver &solver_;

    boost::filesystem::path path_;

    int iterNo_ = 0, fileWriteFrequency_ = 1;
};


#endif
