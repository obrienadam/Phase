#include <iostream>
#include <memory>

#include "Input.h"
#include "CommandLine.h"
#include "ConstructGrid.h"
#include "TecplotViewer.h"
#include "Poisson.h"

int main(int argc, const char* argv[])
{
    using namespace std;

    Input input;
    CommandLine cl(argc, argv, input);

    input.parseInputFile();

    shared_ptr<FiniteVolumeGrid2D> gridPtr = constructGrid(input);

    TecplotViewer viewer(*gridPtr, input);
    Poisson solver(*gridPtr, input);

    viewer.addFieldToOutput(solver.phi);
    viewer.write(0.);

    Scalar error = solver.solve(0.1);

    cout << "Error: " << error << "\n";

    return 0;
}
