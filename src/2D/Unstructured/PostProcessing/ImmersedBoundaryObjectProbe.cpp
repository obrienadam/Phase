#include "FiniteVolumeGrid2D/BilinearInterpolator.h"

#include "ImmersedBoundaryObjectProbe.h"

ImmersedBoundaryObjectProbe::ImmersedBoundaryObjectProbe(int fileWriteFreq,
                                                         const std::weak_ptr<const ImmersedBoundaryObject> &ibObj,
                                                         const std::weak_ptr<const ScalarFiniteVolumeField> &field,
                                                         const Vector2D &probePos)
        :
        Object(fileWriteFreq),
        ibObj_(ibObj),
        field_(field)
{

    path_ /= "ImmersedBoundaryObjectProbes/" + ibObj_.lock()->name() + "/" + field_.lock()->name();

    if (field_.lock()->grid()->comm().isMainProc())
        createOutputDirectory();

    probePos_ = probePos;

    if (field_.lock())
    {
        if (field_.lock()->grid()->comm().isMainProc())
        {
            std::ofstream fout(getFilename());
            fout << "time,value\n";
            fout.close();
        }
    }
}

void ImmersedBoundaryObjectProbe::compute(Scalar time, bool force)
{
    if (do_update() || force)
    {
        auto field = field_.lock();

        BilinearInterpolator bi(field->grid(), ibObj_.lock()->position() + probePos_);

        if (bi.isValid())
        {
            std::ofstream fout(getFilename(), std::ofstream::app);
            fout << time << "," << bi(*field_.lock()) << "\n";
            fout.close();
        }
    }
}

std::string ImmersedBoundaryObjectProbe::getFilename() const
{
    return (path_ / (ibObj_.lock()->name() + "_" + field_.lock()->name() + "_probe.csv")).string();
}
