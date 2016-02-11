#include "MatplotlibViewer.h"

MatplotlibViewer::MatplotlibViewer()
{
    Py_Initialize();
    matplotlib_ = PyImport_Import(PyString_FromString("matplotlib.pyplot"));
}

MatplotlibViewer::~MatplotlibViewer()
{
    Py_Finalize();
}

void MatplotlibViewer::writeResidual(int iterationNo, Scalar residual)
{
    PyObject *pResidual = PyFloat_FromDouble(residual);
    PyObject *pIterationNo = PyInt_FromLong(long(iterationNo));

    PyList_Append(iterationNos_, pIterationNo);
    PyList_Append(residuals_, pResidual);
}
