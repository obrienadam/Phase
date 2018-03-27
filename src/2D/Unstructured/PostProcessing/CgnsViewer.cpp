#include <sstream>
#include <iomanip>
#include <stdio.h>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <cgnslib.h>


#include "System/Exception.h"

#include "CgnsViewer.h"

CgnsViewer::CgnsViewer(const Input &input, const Solver &solver)
        :
        Viewer(input, solver)
{
    char filename[256];
    sprintf(filename, "solution/Proc%d/", solver.grid()->comm().rank());

    boost::filesystem::create_directories(filename);

    sprintf(filename, "solution/Proc%d/Grid.cgns", solver.grid()->comm().rank());

    gridfile_ = filename;

    int fid, bid, zid, sid;
    cg_open(filename, CG_MODE_WRITE, &fid);
    bid = createBase(fid, filename_);
    zid = createZone(fid, bid, *solver_.grid(), "Cells");
    writeCoords(fid, bid, zid, *solver_.grid());
    writeConnectivity(fid, bid, zid, *solver_.grid());
    writeBoundaryConnectivity(fid, bid, zid, *solver_.grid());
    //writeImmersedBoundaries(fid, solver_);

    //- Write data necessary to reuse partitioned grid
    cg_sol_write(fid, bid, zid, "Info", CGNS_ENUMV(CellCenter), &sid);

    std::vector<int> procNo(solver.grid()->cells().size(), solver.grid()->comm().rank());
    if (!solver.grid()->bufferZones().empty())
        for (int proc = 0; proc < solver.grid()->comm().nProcs(); ++proc)
            for (const Cell &cell: solver.grid()->bufferZones()[proc])
                procNo[cell.id()] = proc;

    int fieldId;
    cg_field_write(fid, bid, zid, sid, CGNS_ENUMV(Integer), "ProcNo", procNo.data(), &fieldId);

    std::vector<int> globalId(solver.grid()->cells().size(), -1);
    std::vector<int> globalIdStart(1, 0);
    CellGroup localCells = solver.grid()->localActiveCells() + solver.grid()->localInactiveCells();
    for (int size: solver.grid()->comm().allGather(localCells.size()))
        globalIdStart.push_back(globalIdStart.back() + size);

    int id = globalIdStart[solver.grid()->comm().rank()];
    for (const Cell &cell: localCells)
        globalId[cell.id()] = id++;

    solver.grid()->sendMessages(globalId);
    cg_field_write(fid, bid, zid, sid, CGNS_ENUMV(Integer), "GlobalID", globalId.data(), &fieldId);
    cg_close(fid);
}

void CgnsViewer::write(Scalar solutionTime)
{
    char filename[256];
    sprintf(filename, "solution/%lf/Proc%d", solutionTime, solver_.grid()->comm().rank());

    boost::filesystem::create_directories(filename);

    sprintf(filename, "solution/%lf/Proc%d/Solution.cgns", solutionTime, solver_.grid()->comm().rank());

    int fid, bid, zid, sid;

    cg_open(filename, CG_MODE_WRITE, &fid);

    bid = createBase(fid, filename_);
    zid = createZone(fid, bid, *solver_.grid(), "Cells");
    linkGrid(fid, bid, zid, solver_.grid()->comm());

    cg_sol_write(fid, bid, zid, "Solution", CGNS_ENUMV(CellCenter), &sid);

    int fieldId;

    for (const FiniteVolumeField<int> &field: integerFields_)
        cg_field_write(fid, bid, zid, sid, CGNS_ENUMV(Integer), field.name().c_str(), field.data(), &fieldId);

    for (const ScalarFiniteVolumeField &field: scalarFields_)
        cg_field_write(fid, bid, zid, sid, CGNS_ENUMV(RealDouble), field.name().c_str(), field.data(), &fieldId);

    for (const VectorFiniteVolumeField &field: vectorFields_)
    {
        std::vector<Scalar> x(field.grid()->nCells()), y(field.grid()->nCells());
        std::transform(field.begin(), field.end(), x.begin(), [](const Vector2D &vec) { return vec.x; });
        std::transform(field.begin(), field.end(), y.begin(), [](const Vector2D &vec) { return vec.y; });

        cg_field_write(fid, bid, zid, sid, CGNS_ENUMV(RealDouble), (field.name() + "X").c_str(), x.data(), &fieldId);
        cg_field_write(fid, bid, zid, sid, CGNS_ENUMV(RealDouble), (field.name() + "Y").c_str(), y.data(), &fieldId);
    }

    cg_close(fid);
}

int CgnsViewer::createBase(int fid, const std::string &name)
{
    int id;

    if (name.size() > 32)
        throw Exception("CgnsViewer", "createBase",
                        "case name \"" + name + "\" is too long, names can only contain a max of 32 characters.");

    cg_base_write(fid, name.c_str(), 2, 2, &id);
    //cg_simulation_type_write(fid, id, TimeAccurate);

    return id;
}

int CgnsViewer::createZone(int fid, int bid, const FiniteVolumeGrid2D &grid, const std::string &name)
{
    cgsize_t sizes[3] = {(cgsize_t) grid.nNodes(), (cgsize_t) grid.nCells(), 0};
    int id;
    cg_zone_write(fid, bid, name.c_str(), sizes, CGNS_ENUMV(Unstructured), &id);

    return id;
}

void CgnsViewer::writeCoords(int fid, int bid, int zid, const FiniteVolumeGrid2D &grid)
{
    //- Write the grid info
    std::vector<Scalar> coordsX = grid.xCoords(), coordsY = grid.yCoords();

    int xid;
    cg_coord_write(fid, bid, zid, CGNS_ENUMV(RealDouble), "CoordinateX", coordsX.data(), &xid);
    cg_coord_write(fid, bid, zid, CGNS_ENUMV(RealDouble), "CoordinateY", coordsY.data(), &xid);
}

int CgnsViewer::writeConnectivity(int fid, int bid, int zid, const FiniteVolumeGrid2D &grid)
{
    //- Write connectivity
    std::vector<cgsize_t> connectivity;
    connectivity.reserve(5 * grid.nCells());

    for (const Cell &cell: grid.cells())
    {
        switch (cell.nodes().size())
        {
            case 3:
                connectivity.push_back(CGNS_ENUMV(TRI_3));
                break;
            case 4:
                connectivity.push_back(CGNS_ENUMV(QUAD_4));
                break;
        }

        for (const Node &node: cell.nodes())
            connectivity.push_back(node.id() + 1);
    }

    int id;
    cg_section_write(fid, bid, zid, "GridElements", CGNS_ENUMV(MIXED), 1, grid.cells().size(), 0, connectivity.data(),
                     &id);

    return id;
}

void CgnsViewer::writeBoundaryConnectivity(int fid, int bid, int zid, const FiniteVolumeGrid2D &grid)
{
    //- Now write the boundary mesh elements
    cgsize_t start = grid.nCells() + 1;
    for (const Patch &patch: grid.patches())
    {
        cgsize_t end = start + patch.size() - 1;
        std::vector<cgsize_t> connectivity;

        std::vector<cgsize_t> elemIds;
        cgsize_t elemId = start;
        for (const Face &face: patch)
        {
            connectivity.push_back(face.lNode().id() + 1);
            connectivity.push_back(face.rNode().id() + 1);
            elemIds.push_back(elemId++);
        }

        int secId;
        cg_section_write(fid, bid, zid, (patch.name() + "Elements").c_str(), CGNS_ENUMV(BAR_2), start, end, 0,
                         connectivity.data(), &secId);

        int bcId;
        cg_boco_write(fid, bid, zid, patch.name().c_str(), CGNS_ENUMV(BCGeneral), CGNS_ENUMV(PointList), elemIds.size(),
                      elemIds.data(), &bcId);
        cg_boco_gridlocation_write(fid, bid, zid, bcId, CGNS_ENUMV(EdgeCenter));

        start = end + 1;
    }
}

void CgnsViewer::linkGrid(int fid, int bid, int zid, const Communicator &comm)
{
    char filename[256];

    sprintf(filename, "../../Proc%d/Grid.cgns", comm.rank());

    cg_goto(fid, bid, "Zone_t", zid, "end");
    cg_link_write("GridCoordinates", filename, ("/" + filename_ + "/Cells/GridCoordinates").c_str());
    cg_link_write("GridElements", filename, ("/" + filename_ + "/Cells/GridElements").c_str());

    if (!solver_.grid()->patches().empty())
        cg_link_write("ZoneBC", filename, ("/" + filename_ + "/Cells/ZoneBC").c_str());

    for (const Patch &patch: solver_.grid()->patches())
    {
        std::string name(patch.name());
        cg_link_write((name + "Elements").c_str(), filename, ("/" + filename_ + "/Cells/" + name + "Elements").c_str());
    }
}
