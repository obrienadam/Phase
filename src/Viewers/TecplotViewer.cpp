#include "TecplotViewer.h"

TecplotViewer::TecplotViewer(const Input &input, const Communicator &comm, const Solver &solver)
        :
        Viewer(input, comm, solver)
{
    if (comm.isMainProc())
    {
        for (const ImmersedBoundaryObject &ibObj: solver_.ibObjs())
            outFiles_.push_back(std::ofstream((ibObj.name() + ".dat").c_str()));

        auto fItr = outFiles_.begin();
        for (const ImmersedBoundaryObject &ibObj: solver_.ibObjs())
            *fItr++ << "Title = \"" << ibObj.name() << "\"\n";
    }
}

TecplotViewer::~TecplotViewer()
{
    close();
}

void TecplotViewer::close()
{
    for (std::ofstream &fout: outFiles_)
        fout.close();
}

void TecplotViewer::write(Scalar solutionTime, const Communicator &comm)
{
    if (solver_.comm().isMainProc())
    {
        auto fItr = outFiles_.begin();

        for (const ImmersedBoundaryObject &ibObj: solver_.ibObjs())
        {
            const Shape2D *shape = &ibObj.shape();

            //- Add a geometry record for this zone
            switch (ibObj.shape().type())
            {
                case Shape2D::CIRCLE:
                {
                    const Circle &circ = *static_cast<const Circle*>(shape);

                    *fItr << "Geometry x=" << circ.centroid().x << " y=" << circ.centroid().y
                          << " T=CIRCLE C=BLACK LT=0.4 FC=CUST2 CS=GRID EP=300 ZN=" << zoneNo_ << "\n"
                          << circ.radius() << std::endl;
                }
                    break;

                case Shape2D::BOX:
                {
                    const Box& box = *static_cast<const Box*>(shape);

                    *fItr << "Geometry x=" << box.lower().x << " y=" << box.lower().y
                          << " T=RECTANGLE C=BLACK LT=0.4 FC=CUST2 CS=GRID ZN=" << zoneNo_ << "\n"
                          << box.upper().x - box.lower().x << " " << box.upper().y - box.lower().y << std::endl;
                }
                    break;
                default:
                {
                    auto pgn = shape->polygonize();

                    *fItr << "Geometry x=0 y=0 T=LINE C=BLACK LT=0.4 FC=CUST2 CS=GRID ZN=" << zoneNo_ << "\n"
                          << "1\n"
                          << pgn.vertices().size() << "\n";

                    for (const Point2D& pt: pgn.vertices())
                        *fItr << pt.x << " " << pt.y << "\n";

                    fItr->flush();
                }

                    break;
            }

            *fItr++;
        }
    }

    zoneNo_++;
}
