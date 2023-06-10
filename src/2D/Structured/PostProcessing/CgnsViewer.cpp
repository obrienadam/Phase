#include <boost/algorithm/string.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

#include "System/CgnsFile.h"

#include "CgnsViewer.h"

CgnsViewer::CgnsViewer(const Input &input, const Solver &solver)
    : _solver(solver) {
  boost::filesystem::path path =
      "solution/Proc" + std::to_string(_solver.grid()->comm().rank());

  boost::filesystem::create_directories(path);

  std::string gridfile = (path / "Grid.cgns").string();

  CgnsFile file(gridfile, CgnsFile::WRITE);

  int bid = file.createBase("Grid", 2, 2);
  int zid = file.createStructuredZone(
      bid, "Zone", _solver.grid()->nNodesI(), _solver.grid()->nNodesJ(),
      _solver.grid()->nCellsI(), _solver.grid()->nCellsJ());

  file.writeCoordinates(bid, zid, _solver.grid()->nodes());

  int sid = file.writeSolution(bid, zid, "Info");

  file.writeField(bid, zid, sid, "ProcNo", _solver.grid()->ownership());

  std::vector<Label> ids(_solver.grid()->cells().size());

  std::transform(_solver.grid()->cells().begin(), _solver.grid()->cells().end(),
                 ids.begin(), [](const Cell &c) { return c.gid(); });

  file.writeField(bid, zid, sid, "GlobalID", ids);

  std::transform(_solver.grid()->cells().begin(), _solver.grid()->cells().end(),
                 ids.begin(), [](const Cell &c) { return c.lid(); });

  file.writeField(bid, zid, sid, "LocalID", ids);

  file.close();

  std::string scalarFields =
      input.caseInput().get<std::string>("Viewer.scalarFields", "");
  std::string vectorFields =
      input.caseInput().get<std::string>("Viewer.vectorFields", "");

  // split(integerFields_, integerFields, is_any_of(", "), token_compress_on);
  boost::algorithm::split(_scalarFields, scalarFields,
                          boost::algorithm::is_any_of(", "),
                          boost::algorithm::token_compress_on);
  boost::algorithm::split(_vectorFields, vectorFields,
                          boost::algorithm::is_any_of(", "),
                          boost::algorithm::token_compress_on);
}

void CgnsViewer::write(Scalar time) {}
