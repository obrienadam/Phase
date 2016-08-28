#include <sstream>
#include <iomanip>

#include <boost/algorithm/string.hpp>
#include <cgnslib.h>

#include "Viewer.h"
#include "Exception.h"

Viewer::Viewer(const Solver &solver, const Input &input)
    :
      solver_(solver),
      caseName_(input.caseInput().get<std::string>("CaseName"))
{
    using namespace std;
    using namespace boost;

    outputFilename_ = input.outputPath + "/" + caseName_ + ".cgns";

    string vectorFields = input.caseInput().get<string>("Viewer.vectorFields");
    string scalarFields = input.caseInput().get<string>("Viewer.scalarFields");

    vector<string> vectorFieldNames, scalarFieldNames;
    split(vectorFieldNames, vectorFields, is_any_of(", "), token_compress_on);
    split(scalarFieldNames, scalarFields, is_any_of(", "), token_compress_on);

    for(const auto& field: solver_.scalarFields())
    {
        if(std::find(scalarFieldNames.begin(), scalarFieldNames.end(), field.first) != scalarFieldNames.end())
            scalarFields_.push_back(Ref<const ScalarFiniteVolumeField>(field.second));
    }

    for(const auto& field: solver_.vectorFields())
    {
        if(std::find(vectorFieldNames.begin(), vectorFieldNames.end(), field.first) != vectorFieldNames.end())
            vectorFields_.push_back(Ref<const VectorFiniteVolumeField>(field.second));
    }

    cg_open(outputFilename_.c_str(), CG_MODE_WRITE, &fileId_);
    cg_base_write(fileId_, caseName_.c_str(), 2, 2, &baseId_);
    cg_simulation_type_write(fileId_, baseId_, TimeAccurate);

    std::vector<int> sizes;
    sizes.push_back(solver.grid().nNodes());
    sizes.push_back(solver.grid().nCells());
    sizes.push_back(0);

    cg_zone_write(fileId_, baseId_, "Zone 1", sizes.data(), Unstructured, &zoneId_);

    //- Write the grid info
    std::vector<Scalar> coordsX, coordsY;
    coordsX.reserve(solver.grid().nNodes());
    coordsY.reserve(solver.grid().nNodes());

    for(const Node& node: solver.grid().nodes())
    {
        coordsX.push_back(node.x);
        coordsY.push_back(node.y);
    }

    int xid;
    cg_coord_write(fileId_, baseId_, zoneId_, RealDouble, "CoordinateX", coordsX.data(), &xid);
    cg_coord_write(fileId_, baseId_, zoneId_, RealDouble, "CoordinateY", coordsY.data(), &xid);

    //- Write connectivity
    std::vector<cgsize_t> connectivity;
    connectivity.reserve(5*solver.grid().nNodes());

    for(const Cell& cell: solver.grid().cells())
    {
        switch(cell.nodes().size())
        {
        case 3:
            connectivity.push_back(TRI_3);
            break;
        case 4:
            connectivity.push_back(QUAD_4);
            break;
        }

        for(const Node& node: cell.nodes())
            connectivity.push_back(node.id() + 1);
    }

    int secId;
    cg_section_write(fileId_, baseId_, zoneId_, "GridElements", MIXED, 1, solver.grid().cells().size(), 0, connectivity.data(), &secId);

    //- Now write the boundary mesh elements
    cgsize_t start = solver.grid().nCells() + 1;
    for(const auto& patchEntry: solver.grid().patches())
    {
        const Patch& patch = patchEntry.second;

        cgsize_t end = start + patch.faces().size() - 1;
        connectivity.clear();

        std::vector<cgsize_t> elemIds;
        cgsize_t elemId = start;
        for(const Face &face: patch.faces())
        {
            connectivity.push_back(face.lNode().id() + 1);
            connectivity.push_back(face.rNode().id() + 1);
            elemIds.push_back(elemId++);
        }

        cg_section_write(fileId_, baseId_, zoneId_, (patch.name + "Elements").c_str(), BAR_2, start, end, 0, connectivity.data(), &secId);

        int bcId;
        cg_boco_write(fileId_, baseId_, zoneId_, patch.name.c_str(), BCGeneral, PointList, elemIds.size(), elemIds.data(), &bcId);
        cg_boco_gridlocation_write(fileId_, baseId_, zoneId_, bcId, EdgeCenter);

        start = end + 1;
    }

    cg_close(fileId_);
}

void Viewer::write(Scalar solutionTime)
{
    const std::string solutionName = "FlowSolution" + std::to_string(flowSolutionPointers_.size() + 1);

    //- Update time step info
    timeValues_.push_back(solutionTime);
    flowSolutionPointers_.push_back(solutionName);

    //- Open the file
    cg_open(outputFilename_.c_str(), CG_MODE_MODIFY, &fileId_);

    //- Rewrite the iteration info (necessary unfortunately)
    cg_biter_write(fileId_, baseId_, "TimeIterValues", timeValues_.size());
    cg_goto(fileId_, baseId_, "BaseIterativeData_t", 1, "end");

    cgsize_t dim[2];
    dim[0] = timeValues_.size();

    cg_array_write("TimeValues", RealDouble, 1, dim, timeValues_.data());

    //- Write the solution data for the current time step
    int solutionId;
    cg_sol_write(fileId_, baseId_, zoneId_, solutionName.c_str(), CellCenter, &solutionId);

    int fieldId;
    for(const ScalarFiniteVolumeField& field: scalarFields_)
        cg_field_write(fileId_, baseId_, zoneId_, solutionId, RealDouble, field.name.c_str(), field.data(), &fieldId);

    for(const VectorFiniteVolumeField& field: vectorFields_)
    {
        std::vector<Scalar> xComps, yComps;

        xComps.reserve(field.size());
        yComps.reserve(field.size());

        for(const Vector2D& vec: field)
        {
            xComps.push_back(vec.x);
            yComps.push_back(vec.y);
        }

        cg_field_write(fileId_, baseId_, zoneId_, solutionId, RealDouble, (field.name + "X").c_str(), xComps.data(), &fieldId);
        cg_field_write(fileId_, baseId_, zoneId_, solutionId, RealDouble, (field.name + "Y").c_str(), yComps.data(), &fieldId);
    }

    //- Update the flow solution pointers

    cg_ziter_write(fileId_, baseId_, zoneId_, "ZoneIterativeData");
    cg_goto(fileId_, baseId_, "Zone_t", zoneId_, "ZoneIterativeData_t", 1, "end");

    std::ostringstream sout;

    for(const std::string& str: flowSolutionPointers_)
    {
        sout.width(32); // Interestingly enough, this value must be 32
        sout << std::left << std::setfill(' ') << str;
    }

    dim[0] = 32;
    dim[1] = flowSolutionPointers_.size();

    cg_array_write("FlowSolutionPointers", Character, 2, dim, sout.str().c_str());

    cg_close(fileId_);
}
