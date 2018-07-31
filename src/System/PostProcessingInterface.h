#ifndef PHASE_POST_PROCESSING_INTERFACE_H
#define PHASE_POST_PROCESSING_INTERFACE_H

#include <memory>

#include <boost/filesystem.hpp>

#include "Types/Types.h"

class PostProcessingInterface
{

public:

    class Object
    {
    public:

        Object(int fileWriteFreq);

        virtual void compute(Scalar timeStep) = 0;

        virtual bool do_update();

    protected:

        void createOutputDirectory() const;

        boost::filesystem::path path_;

        int iter_ = 0, fileWriteFreq_ = 1;
    };

    virtual void compute(Scalar time);

protected:

    boost::filesystem::path path_;

    std::vector<std::shared_ptr<Object>> objs_;
};

#endif
