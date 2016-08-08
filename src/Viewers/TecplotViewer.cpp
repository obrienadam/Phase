#include "TecplotViewer.h"
#include "Exception.h"

TecplotViewer::TecplotViewer(const Solver &solver, const Input &input)
    :
      Viewer(solver, input)
{
    createTecplotHeader();
    nZones_ = 0;

    for(const auto& entry: solver.geometries())
        geometryRecords_.insert(std::make_pair(entry.first, std::cref(entry.second)));
}

void TecplotViewer::createTecplotHeader()
{   
    openFile();

    fout_ << "Title = \"" << caseName_ << "\"\n"
          << "Variables = \"x\", \"y\"";

    for(const ScalarFiniteVolumeField& field: scalarFields_)
    {
        fout_ << ", \"" << field.name << "\"";
    }

    for(const VectorFiniteVolumeField& field: vectorFields_)
    {
        fout_ << ", \"" << field.name << "_x\", \"" << field.name << "_y\"";
    }

    fout_ << "\n";
}

void TecplotViewer::write(Scalar solutionTime)
{
    fout_ << "Zone T = \"" << caseName_ << "_" << solutionTime << "s\"\n"
          << "N = " << solver_.grid().nodes().size() << ", E = " << solver_.grid().cells().size() << ", ZoneType = FeQuadrilateral, Datapacking = Block\n"
          << "Varlocation=(" << varLocation() << "=CellCentered)\n"
          << "StrandId = 1, SolutionTime = " << solutionTime << "\n";

    for(const Node& node: solver_.grid().nodes())
        fout_ << node.x << "\n";

    for(const Node& node: solver_.grid().nodes())
        fout_ << node.y << "\n";

    for(const ScalarFiniteVolumeField& field: scalarFields_)
        for(Scalar val: field)
            fout_ << val << "\n";

    for(const VectorFiniteVolumeField& field: vectorFields_)
    {
        for(const Vector2D& vec: field)
            fout_ << vec.x << "\n";

        for(const Vector2D& vec: field)
            fout_ << vec.y << "\n";
    }

    for(const Cell& cell: solver_.grid().cells())
    {
        for(const Node& node: cell.nodes())
            fout_ << node.id() + 1 << " ";
        fout_ << "\n";
    }

    ++nZones_;
    //- Create geometry records
    for(const auto& entry: geometryRecords_)
    {
        const std::vector<Polygon>& pgns = entry.second;

        for(const Polygon &pgn: pgns)
        {
            if(pgn.vertices().size() == 0)
                continue;

            Point2D origin = *pgn.begin();

            fout_ << "Geometry x=" << origin.x << ", y=" << origin.y << ", t=line, fc=blue, cs=grid, zn=" << nZones_ << "\n"
                  << "1\n"
                  << pgn.vertices().size() << "\n";

            for(const Point2D& vtx: pgn)
            {
                Vector2D rVec = vtx - origin;
                fout_ << rVec.x << " " << rVec.y << "\n";
            }
        }
    }
}

std::string TecplotViewer::varLocation()
{
    size_t nCellCenteredVariables = scalarFields_.size() + 2*vectorFields_.size();
    return nCellCenteredVariables == 1 ? "[3]" : nCellCenteredVariables == 2 ? "[3,4]" : "[3-" + std::to_string(nCellCenteredVariables + 2) + "]";
}
