#ifndef MATPLOTLIB_VIEWER_H
#define MATPLOTLIB_VIEWER_H

#include <vector>

#include <python2.7/Python.h>

#include "Types.h"

class MatplotlibViewer
{
public:
    MatplotlibViewer();
    ~MatplotlibViewer();

    void writeResidual(int iterationNo, Scalar residual);

protected:

    PyObject* matplotlib_, *iterationNos_, *residuals_;
};

#endif
