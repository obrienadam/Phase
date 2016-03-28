#include "TecplotViewer.h"
#include "Exception.h"

TecplotViewer::TecplotViewer(const Solver &solver, const Input &input)
    :
      Viewer(solver, input)
{
    createTecplotHeader();
}

void TecplotViewer::createTecplotHeader()
{   
    openFile();

    fout_ << "Title = \"" << caseName_ << "\"\n"
          << "Variables = \"x\", \"y\"";

    for(const auto& field: grid_.scalarFields())
    {
        fout_ << ", \"" << field.first << "\"";
        scalarFields_.push_back(Ref<const ScalarFiniteVolumeField>(field.second));
    }

    for(const auto& field: grid_.vectorFields())
    {
        fout_ << ", \"" << field.first << "_x\", \"" << field.first << "_y\"";
    }

    fout_ << "\n";
}

void TecplotViewer::write(Scalar solutionTime)
{
    fout_ << "Zone T = \"" << caseName_ << "_" << solutionTime << "s\"\n"
          << "N = " << grid_.nodes.size() << ", E = " << grid_.cells.size() << ", ZoneType = FeQuadrilateral, Datapacking = Block\n"
          << "Varlocation=(" << varLocation() << "=CellCentered)\n";

    for(const Node& node: grid_.nodes)
        fout_ << node.x << "\n";

    for(const Node& node: grid_.nodes)
        fout_ << node.y << "\n";

    for(const ScalarFiniteVolumeField& field: scalarFields_)
        for(Scalar val: field)
            fout_ << val << "\n";

    for(const Cell& cell: grid_.cells)
    {
        const auto &nodeIds = cell.nodeIds();
        fout_ << nodeIds[0] + 1 << " " << nodeIds[1] + 1 << " " << nodeIds[2] + 1 << " " << nodeIds[3] + 1 << "\n";
    }
}

std::string TecplotViewer::varLocation()
{
    size_t nCellCenteredVariables = scalarFields_.size() + 2*vectorFields_.size();
    return nCellCenteredVariables == 1 ? "[3]" : nCellCenteredVariables == 2 ? "[3,4]" : "[3-" + std::to_string(nCellCenteredVariables + 2) + "]";
}
