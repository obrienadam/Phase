#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#include "Solvers/Solver.h"
#include "CompactCgnsViewer.h"

CompactCgnsViewer::CompactCgnsViewer(const CommandLine &cl, const Input &input, const Solver &solver)
    :
      Viewer(cl, input, solver),
      solnNo_(0)
{
    boost::filesystem::path path = "./solution";

    if(solver.comm().isMainProc())
        boost::filesystem::create_directory(path);

    solver.comm().barrier();

    if(!isRestart_)
    {
        filename_ = (path / ("proc" + std::to_string(solver.comm().rank()) + ".cgns")).string();
        CgnsFile file(filename_, CgnsFile::WRITE);

        bid_ = file.createBase(input.caseInput().get<std::string>("CaseName"), 2, 2);
        zid_ = file.createUnstructuredZone(bid_, "Cells", solver.grid()->nNodes(), solver.grid()->nCells());

        file.writeCoordinates(bid_, zid_, solver.grid()->coords());

        auto cptr = solver.grid()->eptr();
        auto cind = solver.grid()->eind();

        std::transform(cind.begin(), cind.end(), cind.begin(), [](Label id)
        { return id + 1; });

        file.writeMixedElementSection(bid_, zid_, "Elements", 1, solver.grid()->nCells(), cptr, cind);

        //- Now write the boundary mesh elements
        //    size_t start = solver.grid()->nCells() + 1;
        //    for (const FaceGroup &patch: solver.grid()->patches())
        //    {
        //        size_t end = start + patch.size() - 1;

        //        std::vector<int> elems;

        //        for (const Face &face: patch)
        //            elems.insert(elems.end(), {(int)face.lNode().id() + 1, (int)face.rNode().id() + 1});

        //        int sid = file.writeBarElementSection(bid_, zid_, (patch.name() + "Elements"), start, end, elems);
        //        int bcid = file.writeBoCo(bid_, zid_, patch.name(), start, end);

        //        start = end + 1;
        //    }

        //- Domain info
        int sid = file.writeSolution(bid_, zid_, "Info");
        file.writeField(bid_, zid_, sid, "ProcNo", solver.grid()->cellOwnership());
        file.writeField(bid_, zid_, sid, "GlobalID", solver.grid()->globalIds());
        file.close();
    }
}

void CompactCgnsViewer::write(Scalar time)
{
    CgnsFile file(filename_, CgnsFile::MODIFY);

    int sid = file.writeSolution(bid_, zid_, "FlowSolution" + std::to_string(++solnNo_));
    file.writeDescriptorNode(bid_, zid_, sid, "SolutionTime", std::to_string(time));

    for (const std::string &fieldname: integerFields_)
    {
        auto field = solver_.integerField(fieldname);
        if (field)
            file.writeField(bid_, zid_, sid, field->name(), *field);
    }

    for (const std::string &fieldname: scalarFields_)
    {
        auto field = solver_.scalarField(fieldname);
        if (field)
            file.writeField(bid_, zid_, sid, field->name(), *field);
    }

    for (const std::string &fieldname: vectorFields_)
    {
        auto field = solver_.vectorField(fieldname);
        if (field)
            file.writeField(bid_, zid_, sid, field->name(), *field);
    }

    file.close();
}
