#ifndef TECPLOT_VIEWER_H
#define TECPLOT_VIEWER_H

#include <fstream>
#include <vector>

#include "Viewer.h"

class TecplotViewer : public Viewer
{
public:

    TecplotViewer(const Solver& solver, const Input& input);

    void write(Scalar solutionTime);

private:

    std::string varLocation();
    void createTecplotHeader();

    std::map<std::string, Ref<const std::vector<Polygon> > > geometryRecords_;
    int nZones_;

};

#endif
