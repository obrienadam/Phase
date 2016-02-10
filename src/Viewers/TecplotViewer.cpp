#include "TecplotViewer.h"
#include "Exception.h"

TecplotViewer::TecplotViewer(const FiniteVolumeGrid2D &grid, const Input &input)
    :
      grid_(grid),
      caseName_(input.get<std::string>("CaseName"))
{
    outputFilename_ = input.outputPath + "/" + caseName_ + ".dat";
}

TecplotViewer::~TecplotViewer()
{
    if(fout_.is_open())
        fout_.close();
}

void TecplotViewer::addFieldToOutput(const ScalarFiniteVolumeField &field)
{
    if(fout_.is_open())
        throw Exception("TecplotViewer", "addFieldToOutput", "cannot add fields to output after first write.");

    scalarFields_.push_back(&field);
}

void TecplotViewer::addFieldToOutput(const VectorFiniteVolumeField &field)
{
    if(fout_.is_open())
        throw Exception("TecplotViewer", "addFieldToOutput", "cannot add fields to output after first write.");

    vectorFields_.push_back(&field);
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

void TecplotViewer::writeTec360(Scalar solutionTime)
{
    bool writeMesh = false;

    if(!fout_.is_open())
    {
        writeMesh = true;
        fout_.open(outputFilename_.c_str());
        createTecplotHeader();
    }

    fout_ << "ZONE T = \"" << caseName_ << "_time:" << solutionTime << "s\", I = " << grid_.nNodesI() << ", J = " << grid_.nNodesJ() << ", F = BLOCK\n"
          << "STRANDID = 1, SOLUTIONTIME = " << solutionTime << "\n";

    uint nOutputVariables = scalarFields_.size() + 2*vectorFields_.size();

    switch (nOutputVariables)
    {
    case 0:
        break;
    case 1:
        fout_ << "VARLOCATION = ([3] = CELLCENTERED)\n";
        break;
    default:
        fout_ << "VARLOCATION = ([3-" << nOutputVariables + 2 << "] = CELLCENTERED)\n";
    }

    if(writeMesh)
    {
        for(int component = 0; component < 2; ++component)
            for(int j = 0; j < grid_.nNodesJ(); ++j)
            {
                for(int i = 0; i < grid_.nNodesI(); ++i)
                {
                    fout_ << grid_.cornerNode(i, j)(component) << " ";
                }
                fout_ << "\n";
            }
    }
    else
    {
        fout_ << "VARSHARELIST = ([1-2] = 1)\n";
    }

    for(const ScalarFiniteVolumeField* scalarFieldPtr: scalarFields_)
    {
        const ScalarFiniteVolumeField& scalarField = *scalarFieldPtr;

        for(int j = 0; j < scalarField.sizeJ(); ++j)
        {
            for(int i = 0; i < scalarField.sizeI(); ++i)
            {
                fout_ << scalarField(i, j) << " ";
            }
            fout_ << "\n";
        }
    }

    for(const VectorFiniteVolumeField* vectorFieldPtr: vectorFields_)
    {
        const VectorFiniteVolumeField& vectorField = *vectorFieldPtr;

        for(int component = 0; component < 2; ++component)
            for(int j = 0; j < vectorField.sizeJ(); ++j)
            {
                for(int i = 0; i < vectorField.sizeI(); ++i)
                {
                    fout_ << vectorField(i, j)(component) << " ";
                }
                fout_ << "\n";
            }
    }
}
