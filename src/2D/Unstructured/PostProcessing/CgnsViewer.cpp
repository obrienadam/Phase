#include <sstream>
#include <iomanip>
#include <stdio.h>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

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

    auto cptr = solver.grid()->eptr();
    auto cind = solver.grid()->eind();

    std::transform(cind.begin(), cind.end(), cind.begin(), [](Label id)
    { return id + 1; });

    file.writeMixedElementSection(bid, zid, "Cells", 1, solver.grid()->nCells(), cptr, cind);

    //- Now write the boundary mesh elements
    size_t start = solver.grid()->nCells() + 1;
    for (const FaceGroup &patch: solver.grid()->patches())
    {
        size_t end = start + patch.size() - 1;

        std::vector<int> elems;

        for (const Face &face: patch)
            elems.insert(elems.end(), {(int)face.lNode().id() + 1, (int)face.rNode().id() + 1});

        int sid = file.writeBarElementSection(bid, zid, (patch.name() + "Elements"), start, end, elems);
        int bcid = file.writeBoCo(bid, zid, patch.name(), start, end);

        start = end + 1;
    }

    int sid = file.writeSolution(bid, zid, "Info");

    file.writeField(bid, zid, sid, "ProcNo", solver.grid()->cellOwnership());
    file.writeField(bid, zid, sid, "GlobalID", solver.grid()->globalIds());

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

    for (const std::string &fieldname: integerFields_)
    {
        auto field = solver_.integerField(fieldname);
        if (field)
            file.writeField(bid, zid, sid, field->name(), *field);
    }

    for (const std::string &fieldname: scalarFields_)
    {
        auto field = solver_.scalarField(fieldname);
        if (field)
            file.writeField(bid, zid, sid, field->name(), *field);
    }

    for (const std::string &fieldname: vectorFields_)
    {
        auto field = solver_.vectorField(fieldname);
        if (field)
            file.writeField(bid, zid, sid, field->name(), *field);
    }

    path = boost::filesystem::path("../../../") / gridfile_;

    file.linkNode(bid, zid, "GridCoordinates", path.c_str(), "/Grid/Zone/GridCoordinates");
    file.linkNode(bid, zid, "Cells", path.c_str(), "/Grid/Zone/Cells");

    if(!solver_.grid()->patches().empty())
        file.linkNode(bid, zid, "ZoneBC", path.c_str(), "/Grid/Zone/ZoneBC");

    for (const FaceGroup &patch: solver_.grid()->patches())
    {
        file.linkNode(bid, zid, (patch.name() + "Elements").c_str(),
                      path.c_str(),
                      ("/Grid/Zone/" + patch.name() + "Elements").c_str());
    }

    file.linkNode(bid, zid, sid, "GlobalID", path.c_str(), "/Grid/Zone/Info/GlobalID");
    file.linkNode(bid, zid, sid, "ProcNo", path.c_str(), "/Grid/Zone/Info/ProcNo");
    file.close();
}
