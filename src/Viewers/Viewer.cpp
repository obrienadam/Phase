#include <sstream>
#include <iomanip>

#include <boost/algorithm/string.hpp>
#include <cgnslib.h>

#include "Viewer.h"
#include "Exception.h"

Viewer::Viewer(const Input &input, const Solver &solver)
    :
      solver_(solver)
{
    init(input, solver, input.outputPath + "/"
         + input.caseInput().get<std::string>("CaseName")
         + ".cgns");
}

Viewer::Viewer(const Input &input, const Communicator &comm, const Solver &solver)
    :
      solver_(solver)
{
    init(input, comm, solver, input.outputPath + "/"
         + input.caseInput().get<std::string>("CaseName")
         + ".cgns");
}

void Viewer::write(Scalar solutionTime)
{
    //- Open the file
    cg_open(filename_.c_str(), CG_MODE_MODIFY, &fileId_);

    //- Create a new solution node, get the id
    updateBaseIterativeData(solutionTime);
    int solutionId = updateFlowSolutionPointers();

    int fieldId;
    for(const ScalarFiniteVolumeField& field: scalarFields_)
        cg_field_write(fileId_, baseId_, zoneId_, solutionId, RealDouble, field.name().c_str(), field.data(), &fieldId);

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

        cg_field_write(fileId_, baseId_, zoneId_, solutionId, RealDouble, (field.name() + "X").c_str(), xComps.data(), &fieldId);
        cg_field_write(fileId_, baseId_, zoneId_, solutionId, RealDouble, (field.name() + "Y").c_str(), yComps.data(), &fieldId);
    }

    cg_close(fileId_);
}

void Viewer::write(Scalar solutionTime, const Communicator &comm)
{
    if(comm.isMainProc())
    {
        cg_open(filename_.c_str(), CG_MODE_MODIFY, &fileId_);
        updateBaseIterativeData(solutionTime);
        cg_close(fileId_);
    }

    comm.barrier();
    for(int proc = 0; proc < comm.nProcs(); ++proc)
    {
        if(proc == comm.rank())
        {
            cg_open(filename_.c_str(), CG_MODE_MODIFY, &fileId_);
            int solutionId = updateFlowSolutionPointers();

            int fieldId;
            for(const ScalarFiniteVolumeField& field: scalarFields_)
                cg_field_write(fileId_, baseId_, zoneId_, solutionId, RealDouble, field.name().c_str(), field.data(), &fieldId);

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

                cg_field_write(fileId_, baseId_, zoneId_, solutionId, RealDouble, (field.name() + "X").c_str(), xComps.data(), &fieldId);
                cg_field_write(fileId_, baseId_, zoneId_, solutionId, RealDouble, (field.name() + "Y").c_str(), yComps.data(), &fieldId);
            }

            cg_close(fileId_);
        }

        comm.barrier();
    }
}

void Viewer::write(const std::vector<VolumeIntegrator> &volumeIntegrators)
{
    cg_open(filename_.c_str(), CG_MODE_MODIFY, &fileId_);

    for(const VolumeIntegrator& vi: volumeIntegrators)
    {
        int nSteps = vi.data().size();

        cg_goto(fileId_, baseId_, "end");
        cg_convergence_write(vi.data().size(), "");
        cg_goto(fileId_, baseId_, "ConvergenceHistory_t", 1, "end");
        cg_array_write((vi.field().name() + "VolumeIntegral").c_str(), RealDouble, 1, &nSteps, vi.data().data());
    }

    cg_close(fileId_);
}

void Viewer::init(const Input &input, const Solver &solver, const std::string &filename)
{
    filename_ = filename;

    initFields(input);

    cg_open(filename_.c_str(), CG_MODE_WRITE, &fileId_);

    baseId_ = createBase(input.caseInput().get<std::string>("CaseName"));
    zoneId_ = createZone(solver.grid());

    writeCoords(solver.grid());
    writeConnectivity(solver.grid());
    writeBoundaryConnectivity(solver.grid());

    if(solver.ib().ibObjs().size() > 0)
        writeImmersedBoundaries(solver);

    cg_close(fileId_);
}

void Viewer::init(const Input &input, const Communicator &comm, const Solver &solver, const std::string &filename)
{
    filename_ = filename;

    initFields(input);

    if(comm.isMainProc())
    {
        cg_open(filename_.c_str(), CG_MODE_WRITE, &fileId_);
        baseId_ = createBase(input.caseInput().get<std::string>("CaseName"));
        cg_close(fileId_);
    }

    baseId_ = comm.broadcast(comm.mainProcNo(), baseId_);
    for(int proc = 0; proc < comm.nProcs(); ++proc)
    {
        if(proc == comm.rank())
        {
            cg_open(filename_.c_str(), CG_MODE_MODIFY, &fileId_);

            zoneId_ = createZone(solver.grid(), "CellsProc" + std::to_string(proc));

            writeCoords(solver.grid());
            writeConnectivity(solver.grid());
            writeBoundaryConnectivity(solver.grid());

            cg_close(fileId_);
        }

        comm.barrier();
    }
}

void Viewer::initFields(const Input &input)
{
    using namespace std;
    using namespace boost;

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
}

int Viewer::createBase(const std::string &name)
{
    int id;
    cg_base_write(fileId_, name.c_str(), 2, 2, &id);
    cg_simulation_type_write(fileId_, id, TimeAccurate);

    return id;
}

int Viewer::createZone(const FiniteVolumeGrid2D &grid, const std::string &name)
{
    cgsize_t sizes[3] = {(cgsize_t)grid.nNodes(), (cgsize_t)grid.nCells(), 0};
    int id;
    cg_zone_write(fileId_, baseId_, name.c_str(), sizes, Unstructured, &id);

    return id;
}

void Viewer::writeCoords(const FiniteVolumeGrid2D &grid)
{
    //- Write the grid info
    std::vector<Scalar> coordsX = grid.xCoords(), coordsY = grid.yCoords();

    int xid;
    cg_coord_write(fileId_, baseId_, zoneId_, RealDouble, "CoordinateX", coordsX.data(), &xid);
    cg_coord_write(fileId_, baseId_, zoneId_, RealDouble, "CoordinateY", coordsY.data(), &xid);
}

int Viewer::writeConnectivity(const FiniteVolumeGrid2D &grid)
{
    //- Write connectivity
    std::vector<cgsize_t> connectivity;
    connectivity.reserve(5*grid.nCells());

    for(const Cell& cell: grid.cells())
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
    cg_section_write(fileId_, baseId_, zoneId_, "GridElements", MIXED, 1, grid.cells().size(), 0, connectivity.data(), &secId);

    return secId;
}

void Viewer::writeBoundaryConnectivity(const FiniteVolumeGrid2D &grid)
{
    //- Now write the boundary mesh elements
    cgsize_t start = grid.nCells() + 1;
    for(const auto& patchEntry: grid.patches())
    {
        const Patch& patch = patchEntry.second;

        cgsize_t end = start + patch.faces().size() - 1;
        std::vector<cgsize_t> connectivity;

        std::vector<cgsize_t> elemIds;
        cgsize_t elemId = start;
        for(const Face &face: patch.faces())
        {
            connectivity.push_back(face.lNode().id() + 1);
            connectivity.push_back(face.rNode().id() + 1);
            elemIds.push_back(elemId++);
        }

        int secId;
        cg_section_write(fileId_, baseId_, zoneId_, (patch.name + "Elements").c_str(), BAR_2, start, end, 0, connectivity.data(), &secId);

        int bcId;
        cg_boco_write(fileId_, baseId_, zoneId_, patch.name.c_str(), BCGeneral, PointList, elemIds.size(), elemIds.data(), &bcId);
        cg_boco_gridlocation_write(fileId_, baseId_, zoneId_, bcId, EdgeCenter);

        start = end + 1;
    }
}

void Viewer::writeImmersedBoundaries(const Solver &solver)
{
    //- Create new zones for the ibs and write them
    int ibBase;
    cg_base_write(fileId_, "IBs", 1, 2, &ibBase);
    for(const ImmersedBoundaryObject& ibObj: solver.ib().ibObjs())
    {
        Polygon pgn = ibObj.shape().polygonize();
        int nVerts = pgn.vertices().size() - 1;
        int sizes[3] = {nVerts, 1, 0};

        int ibZoneId;
        cg_zone_write(fileId_, ibBase, ibObj.name().c_str(), sizes, Unstructured, &ibZoneId);

        std::vector<Scalar> coordsX, coordsY;

        for(const Point2D &vtx: pgn.vertices())
        {
            coordsX.push_back(vtx.x);
            coordsY.push_back(vtx.y);
        }

        int xid;
        cg_coord_write(fileId_, ibBase, ibZoneId, RealDouble, "CoordinateX", coordsX.data(), &xid);
        cg_coord_write(fileId_, ibBase, ibZoneId, RealDouble, "CoordinateY", coordsY.data(), &xid);

        std::vector<cgsize_t> connectivity;

        for(int i = 0; i < nVerts; ++i)
        {
            connectivity.push_back(i + 1);
            connectivity.push_back((i + 1)%nVerts + 1);
        }

        int secId;
        cg_section_write(fileId_, ibBase, ibZoneId, "Edges", BAR_2, 1, nVerts, 0, connectivity.data(), &secId);
    }
}

void Viewer::updateBaseIterativeData(Scalar solutionTime)
{
    timeValues_.push_back(solutionTime);

    //- Rewrite the iteration info (necessary unfortunately)
    cg_biter_write(fileId_, baseId_, "TimeIterValues", timeValues_.size());
    cg_goto(fileId_, baseId_, "BaseIterativeData_t", 1, "end");

    cgsize_t dim[2];
    dim[0] = timeValues_.size();

    cg_array_write("TimeValues", RealDouble, 1, dim, timeValues_.data());
}

int Viewer::updateFlowSolutionPointers()
{
    std::string solutionName = "FlowSolution" + std::to_string(flowSolutionPointers_.size() + 1);

    //- Update time step info and flow solution pointer info
    flowSolutionPointers_.push_back(solutionName);

    //- Rewrite flow solution pointers
    cg_ziter_write(fileId_, baseId_, zoneId_, "ZoneIterativeData");
    cg_goto(fileId_, baseId_, "Zone_t", zoneId_, "ZoneIterativeData_t", 1, "end");

    std::ostringstream sout;

    for(const std::string& str: flowSolutionPointers_)
    {
        sout.width(32); // Interestingly enough, this value must be 32
        sout << std::left << std::setfill(' ') << str;
    }

    cgsize_t dim[2];
    dim[0] = 32;
    dim[1] = flowSolutionPointers_.size();

    cg_array_write("FlowSolutionPointers", Character, 2, dim, sout.str().c_str());

    //- Create a new node which to write
    int solutionId;
    cg_sol_write(fileId_, baseId_, zoneId_, solutionName.c_str(), CellCenter, &solutionId);

    return solutionId;
}
