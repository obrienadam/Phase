#include <sstream>
#include <iomanip>
#include <stdio.h>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#include "System/Exception.h"
#include "System/CgnsFile.h"

#include "CgnsViewer.h"

CgnsViewer::CgnsViewer(const Input &input, const std::weak_ptr<const Solver> &solver)
    :
      _solver(solver)
{
    boost::filesystem::path path = "solution/Proc"
            + std::to_string(_solver.lock()->grid()->comm().rank());

    boost::filesystem::create_directories(path);

    std::string casename = input.caseInput().get<std::string>("CaseName");
    std::string gridfile = (path / "Grid.cgns").string();

    CgnsFile file(gridfile, CgnsFile::WRITE);

    int bid = file.createBase("Grid", 3, 3);
    int zid = file.createStructuredZone(bid,
                                        "Zone",
                                        _solver.lock()->grid()->nNodesI(),
                                        _solver.lock()->grid()->nNodesJ(),
                                        _solver.lock()->grid()->nNodesK(),
                                        _solver.lock()->grid()->nCellsI(),
                                        _solver.lock()->grid()->nCellsJ(),
                                        _solver.lock()->grid()->nCellsK());

    _scalarFieldNames.insert("phi");

    file.writeCoordinates(bid, zid, _solver.lock()->grid()->nodes());
    file.close();
}

void CgnsViewer::write(Scalar time)
{
    boost::filesystem::path path = "solution/"
            + std::to_string(time)
            + "/Proc"
            + std::to_string(_solver.lock()->grid()->comm().rank());

    boost::filesystem::create_directories(path);

    CgnsFile file((path / "Solution.cgns").string(), CgnsFile::WRITE);

    int bid = file.createBase("Solution", 3, 3);
    int zid = file.createStructuredZone(bid,
                                        "Zone",
                                        _solver.lock()->grid()->nNodesI(),
                                        _solver.lock()->grid()->nNodesJ(),
                                        _solver.lock()->grid()->nNodesK(),
                                        _solver.lock()->grid()->nCellsI(),
                                        _solver.lock()->grid()->nCellsJ(),
                                        _solver.lock()->grid()->nCellsK());

    int sid = file.writeSolution(bid, zid, "Solution");

    for(const auto &field: _solver.lock()->scalarFields())
        if(_scalarFieldNames.find(field.first) != _scalarFieldNames.end())
            file.writeField(bid, zid, sid, field.first, field.second->cellData());

    path = boost::filesystem::path("../../") / ("Proc" + std::to_string(_solver.lock()->grid()->comm().rank())) / "Grid.cgns";
    file.linkNode(bid, zid, "GridCoordinates", path.c_str(), "/Grid/Zone/GridCoordinates");

    //    if(!solver_.grid()->patches().empty())
    //        file.linkNode(bid, zid, "ZoneBC", path.c_str(), "/Grid/Zone/ZoneBC");

    //    for (const FaceGroup &patch: solver_.grid()->patches())
    //    {
    //        file.linkNode(bid, zid, (patch.name() + "Elements").c_str(),
    //                      path.c_str(),
    //                      ("/Grid/Zone/" + patch.name() + "Elements").c_str());
    //    }

    //    file.linkNode(bid, zid, sid, "GlobalID", path.c_str(), "/Grid/Zone/Info/GlobalID");
    //    file.linkNode(bid, zid, sid, "ProcNo", path.c_str(), "/Grid/Zone/Info/ProcNo");
    file.close();
}
