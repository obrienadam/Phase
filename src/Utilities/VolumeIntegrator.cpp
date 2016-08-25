#include <boost/algorithm/string.hpp>

#include "VolumeIntegrator.h"
#include "Solver.h"

//- public static functions

std::vector<VolumeIntegrator> VolumeIntegrator::initVolumeIntegrators(const Input &input, const Solver &solver)
{
    std::vector<VolumeIntegrator> volumeIntegrators;

    boost::optional<std::string> volumeIntegratorInput = input.caseInput().get_optional<std::string>("Integrators.VolumeIntegrators");

    if(!volumeIntegratorInput)
        return volumeIntegrators;

    std::string fieldInput = input.caseInput().get<std::string>("Integrators.VolumeIntegrators.fields");
    std::vector<std::string> fieldNames;

    boost::split(fieldNames, fieldInput, boost::is_any_of(", "), boost::token_compress_on);

    std::string cellGroupInput = input.caseInput().get<std::string>("Integrators.VolumeIntegrators.cellGroups");
    std::vector<std::string> cellGroupNames;

    boost::split(cellGroupNames, cellGroupInput, boost::is_any_of(", "), boost::token_compress_on);

    if(fieldNames.size() != cellGroupNames.size())
        throw Exception("VolumeIntegrator", "initVolumeIntegrators", "must specify one cell group for each field name.");

    for(int i = 0; i < fieldNames.size(); ++i)
    {
        auto fieldIt = solver.scalarFields().find(fieldNames[i]);

        if(fieldIt == solver.scalarFields().end())
            continue;

        const ScalarFiniteVolumeField &field = fieldIt->second;
        const CellGroup &cellGroup = solver.grid().cellGroup(cellGroupNames[i]);

        printf("Initializing a volume integrator for field \"%s\" on cell group \"%s\".\n", field.name.c_str(), cellGroup.name().c_str());

        volumeIntegrators.push_back(
                    VolumeIntegrator(field, cellGroup)
                    );
    }

    return volumeIntegrators;
}

//- constructor

VolumeIntegrator::VolumeIntegrator(const ScalarFiniteVolumeField &field, const CellGroup &cellGroup):
    field_(field),
    cellGroup_(cellGroup)
{

}

//- public methods

Scalar VolumeIntegrator::integrate() const
{
    Scalar volInt = 0.;

    for(const Cell &cell: cellGroup_)
        volInt += field_(cell)*cell.volume();

    printf("Volume integral of field \"%s\" over the region \"%s\": %lf\n", field_.name.c_str(), cellGroup_.name().c_str(), volInt);

    return volInt;
}
