#include "ImmersedBoundaryObjectProbe.h"
#include "BilinearInterpolator.h"

ImmersedBoundaryObjectProbe::ImmersedBoundaryObjectProbe(const Solver &solver,
                                                         const std::string &ibObjName,
                                                         const std::string &fieldName,
                                                         const Vector2D &probePos)
        :
        PostProcessingObject(solver)
{

    path_ /= "ImmersedBoundaryObjectProbes/" + ibObjName + "/" + fieldName;

    if (solver.grid()->comm().isMainProc())
        createOutputDirectory();

    ibObj_ = solver_.ib().ibObj(ibObjName);
    field_ = solver_.scalarField(fieldName);
    probePos_ = probePos;

    if (field_.lock())
    {
        if (solver.grid()->comm().isMainProc())
        {
            std::ofstream fout(getFilename());
            fout << "time,value\n";
            fout.close();
        }
    }
}

void ImmersedBoundaryObjectProbe::compute(Scalar time)
{
    if (iterNo_++ % fileWriteFrequency_ == 0)
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