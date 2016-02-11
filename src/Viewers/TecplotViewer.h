#ifndef TECPLOT_VIEWER_H
#define TECPLOT_VIEWER_H

#include <fstream>
#include <vector>

#include "FiniteVolumeGrid2D.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"
#include "Input.h"

class TecplotViewer
{
public:

    TecplotViewer(const FiniteVolumeGrid2D& grid, const Input& input);
    ~TecplotViewer();

    void addFieldToOutput(const ScalarFiniteVolumeField& field);
    void addFieldToOutput(const VectorFiniteVolumeField& field);

    void write(Scalar solutionTime);

private:

    void createTecplotHeader();

    const FiniteVolumeGrid2D& grid_;
    std::vector<const ScalarFiniteVolumeField*> scalarFields_;
    std::vector<const VectorFiniteVolumeField*> vectorFields_;

    std::string caseName_, outputFilename_;
    std::ofstream fout_;

};

#endif
