#include <sstream>
#include <iomanip>
#include <stdio.h>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <cgnslib.h>

#include "System/Exception.h"
#include "System/CgnsFile.h"

#include "CgnsViewer.h"

CgnsViewer::CgnsViewer(const Input &input, const Solver &solver)
        :
        Viewer(input, solver)
{
    boost::filesystem::path path = "solution/Proc" + std::to_string(solver.grid()->comm().rank());
    boost::filesystem::create_directories(path);

    casename_ = input.caseInput().get<std::string>("CaseName");
    gridfile_ = (path / "Grid.cgns").string();

    CgnsFile file(gridfile_, CgnsFile::WRITE);

    int bid = file.createBase("Grid", 2, 2);

    int zid = file.createUnstructuredZone(bid, "Zone", solver.grid()->nNodes(), solver.grid()->nCells());

    file.writeCoordinates(bid, zid, solver.grid()->coords());

    file.writeMixedElementSection(bid, zid, "Cells", 1, solver.grid()->nCells(),
                                  solver.grid()->nodeElementConnectivity());

    //- Now write the boundary mesh elements
    cgsize_t start = solver.grid()->nCells() + 1;
    for (const Patch &patch: solver.grid()->patches())
    {
        cgsize_t end = start + patch.size() - 1;

        int sid = file.writeBarElementSection(bid, zid, (patch.name() + "Elements"), start, end,
                                              patch.faceConnectivity());
        int bcid = file.writeBC(bid, zid, patch.name().c_str(), start, end);

        start = end + 1;
    }

    int sid = file.writeSolution(bid, zid, "Info");

    file.writeField(bid, zid, sid, "ProcNo", *solver.integerField("proc"));

    file.writeField(bid, zid, sid, "GlobalID", *solver.integerField("globalId"));

    file.close();
}

void CgnsViewer::write(Scalar time)
{
    boost::filesystem::path path = "solution/" + std::to_string(time)
                                   + "/Proc" + std::to_string(solver_.grid()->comm().rank());

    boost::filesystem::create_directories(path);

    CgnsFile file((path / "Solution.cgns").string(), CgnsFile::WRITE);

    int bid = file.createBase("Solution", 2, 2);

    int zid = file.createUnstructuredZone(bid, "Zone", solver_.grid()->nNodes(), solver_.grid()->nCells());

    int sid = file.writeSolution(bid, zid, "Solution");

    for (const FiniteVolumeField<int> &field: integerFields_)
        file.writeField(bid, zid, sid, field.name(), field);

    for (const ScalarFiniteVolumeField &field: scalarFields_)
        file.writeField(bid, zid, sid, field.name(), field);

    for (const VectorFiniteVolumeField &field: vectorFields_)
        file.writeField(bid, zid, sid, field.name(), field);

    path = boost::filesystem::path("../../../") / gridfile_;

    file.linkNode(bid, zid, "GridCoordinates", path.c_str(), "/Grid/Zone/GridCoordinates");
    file.linkNode(bid, zid, "Cells", path.c_str(), "/Grid/Zone/Cells");
    file.linkNode(bid, zid, "ZoneBC", path.c_str(), "/Grid/Zone/ZoneBC");

    for (const Patch &patch: solver_.grid()->patches())
    {
        file.linkNode(bid, zid, (patch.name() + "Elements").c_str(),
                      path.c_str(),
                      ("/Grid/Zone/" + patch.name() + "Elements").c_str());
    }

    file.linkNode(bid, zid, sid, "GlobalID", path.c_str(), "/Grid/Zone/Info/GlobalID");
    file.linkNode(bid, zid, sid, "ProcNo", path.c_str(), "/Grid/Zone/Info/ProcNo");
    file.close();

//    char filename[256];
//    sprintf(filename, "solution/%lf/Proc%d", solutionTime, solver_.grid()->comm().rank());
//
//    boost::filesystem::create_directories(filename);
//
//    sprintf(filename, "solution/%lf/Proc%d/Solution.cgns", solutionTime, solver_.grid()->comm().rank());
//
//    int fid, bid, zid, sid;
//
//    cg_open(filename, CG_MODE_WRITE, &fid);
//
//    bid = createBase(fid, filename_);
//    zid = createZone(fid, bid, *solver_.grid(), "Cells");
//    linkGrid(fid, bid, zid, solver_.grid()->comm());
//
//    cg_sol_write(fid, bid, zid, "Solution", CGNS_ENUMV(CellCenter), &sid);
//
//    int fieldId;
//
//    for (const FiniteVolumeField<int> &field: integerFields_)
//        cg_field_write(fid, bid, zid, sid, CGNS_ENUMV(Integer), field.name().c_str(), field.data(), &fieldId);
//
//    for (const ScalarFiniteVolumeField &field: scalarFields_)
//        cg_field_write(fid, bid, zid, sid, CGNS_ENUMV(RealDouble), field.name().c_str(), field.data(), &fieldId);
//
//    for (const VectorFiniteVolumeField &field: vectorFields_)
//    {
//        std::vector<Scalar> x(field.grid()->nCells()), y(field.grid()->nCells());
//        std::transform(field.begin(), field.end(), x.begin(), [](const Vector2D &vec)
//        { return vec.x; });
//        std::transform(field.begin(), field.end(), y.begin(), [](const Vector2D &vec)
//        { return vec.y; });
//
//        cg_field_write(fid, bid, zid, sid, CGNS_ENUMV(RealDouble), (field.name() + "X").c_str(), x.data(), &fieldId);
//        cg_field_write(fid, bid, zid, sid, CGNS_ENUMV(RealDouble), (field.name() + "Y").c_str(), y.data(), &fieldId);
//    }
//
//    cg_close(fid);
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