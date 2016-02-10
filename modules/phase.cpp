#include <iostream>

#include "Input.h"
#include "CommandLine.h"
#include "EquidistantGrid2D.h"
#include "TecplotViewer.h"
#include "ScalarFiniteVolumeField.h"

int main(int argc, const char* argv[])
{
    using namespace std;

    Input input;
    CommandLine cl(argc, argv, input);

    input.parseInputFile();

    EquidistantGrid2D grid(100, 100, 0.1);
    TecplotViewer viewer(grid, input);
    ScalarFiniteVolumeField phi(grid, "phi"), k(grid, "k");
    VectorFiniteVolumeField u(grid, "u"), gradP(grid, "gradP");

    viewer.addFieldToOutput(phi);
    viewer.addFieldToOutput(k);
    viewer.addFieldToOutput(u);
    viewer.addFieldToOutput(gradP);

    viewer.writeTec360(0.);
    viewer.writeTec360(1.);
    viewer.writeTec360(2.);
    viewer.writeTec360(3.);
    viewer.writeTec360(4.);

    return 0;
}
