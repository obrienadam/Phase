#include "TecplotViewer.h"

TecplotViewer::TecplotViewer(const Input &input, const Communicator &comm, const Solver &solver)
    :
      Viewer(input, comm, solver)
{

}

TecplotViewer::~TecplotViewer()
{

}

void TecplotViewer::close()
{


}

void TecplotViewer::write(Scalar solutionTime, const Communicator &comm)
{

}
