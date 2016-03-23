#include "TecplotViewer.h"
#include "Exception.h"

TecplotViewer::TecplotViewer(const Solver &solver, const Input &input)
    :
      Viewer(solver, input)
{

}

void TecplotViewer::createTecplotHeader()
{   
    fout_ << "TITLE = \"" << caseName_ << "\"\n"
          << "VARIABLES = \"x\", \"y\"";

    for(const ScalarFiniteVolumeField* scalarFieldPtr: scalarFields_)
    {
        fout_ << ", \"" << scalarFieldPtr->name << "\"";
    }

    for(const VectorFiniteVolumeField* vectorFieldPtr: vectorFields_)
    {
        fout_ << ", \"" << vectorFieldPtr->name + "_x\", \"" << vectorFieldPtr->name + "_y\"";
    }

    fout_ << "\n";
}

void TecplotViewer::write(Scalar solutionTime)
{
    bool writeMesh = false;

    if(!fout_.is_open())
    {
        writeMesh = true;
        fout_.open(outputFilename_.c_str());
        createTecplotHeader();
    }

//    fout_ << "ZONE T = \"" << caseName_ << "_time:" << solutionTime << "s\", I = " << grid_.nNodesI() << ", J = " << grid_.nNodesJ() << ", F = BLOCK\n"
//          << "STRANDID = 1, SOLUTIONTIME = " << solutionTime << "\n";

//    uint nOutputVariables = scalarFields_.size() + 2*vectorFields_.size();

//    switch (nOutputVariables)
//    {
//    case 0:
//        break;
//    case 1:
//        fout_ << "VARLOCATION = ([3] = CELLCENTERED)\n";
//        break;
//    default:
//        fout_ << "VARLOCATION = ([3-" << nOutputVariables + 2 << "] = CELLCENTERED)\n";
//    }

//    if(writeMesh)
//    {
//        for(int component = 0; component < 2; ++component)
//            for(int j = 0; j < grid_.nNodesJ(); ++j)
//            {
//                for(int i = 0; i < grid_.nNodesI(); ++i)
//                {
//                    fout_ << grid_.cornerNode(i, j)(component) << " ";
//                }
//                fout_ << "\n";
//            }
//    }
    else
    {
        fout_ << "VARSHARELIST = ([1-2] = 1)\n";
    }
}
