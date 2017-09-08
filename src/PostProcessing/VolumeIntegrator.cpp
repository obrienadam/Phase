#include "VolumeIntegrator.h"

VolumeIntegrator::VolumeIntegrator(const Solver &solver, const std::string &fieldName)
        :
        PostProcessingObject(solver),
        field_(solver.scalarField(fieldName))
{
    outputDir_ = outputDir_ / "VolumeIntegrators";

    if (solver.grid().comm().isMainProc())
        createOutputDirectory();

    std::ofstream fout((outputDir_ / (field_.name() + ".dat")).string());
    fout << "time\tresult\n";
    fout.close();
}

void VolumeIntegrator::compute(Scalar time)
{
    if(iterNo_++ % fileWriteFrequency_ == 0)
    {
        Scalar result = 0.;

        for (const Cell &cell: solver_.grid().localActiveCells())
            result += field_(cell) * cell.volume();

        result = solver_.grid().comm().sum(result);

        if (solver_.grid().comm().isMainProc())
        {
            std::ofstream fout((outputDir_ / (field_.name() + ".dat")).string(),
                               std::ofstream::out | std::ofstream::app);
            fout << time << "\t" << result << "\n";
            fout.close();
        }
    }
}