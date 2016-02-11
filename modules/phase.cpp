#include <iostream>

#include "Input.h"
#include "CommandLine.h"
#include "EquidistantGrid2D.h"
#include "TecplotViewer.h"
#include "Simple.h"

int main(int argc, const char* argv[])
{
    using namespace std;

    Input input;
    CommandLine cl(argc, argv, input);
    EquidistantGrid2D grid(100, 100, 0.1);

    input.parseInputFile();

    TecplotViewer viewer(grid, input);
    Simple solver(grid, input);

    viewer.addFieldToOutput(solver.u);
    viewer.addFieldToOutput(solver.p);

    viewer.writeTec360(0.);
    viewer.writeTec360(1.);
    viewer.writeTec360(2.);

    return 0;
}
