#include "IbViewer.h"

IbViewer::IbViewer(const Input &input, const Communicator &comm, const Solver &solver)
    :
      Viewer(input, comm, solver),
      ib_(solver.ib())
{
    if (comm.isMainProc())
    {
        for (const ImmersedBoundaryObject &ibObj: solver_.ibObjs())
        {
            ibFiles_.push_back(std::ofstream((ibObj.name() + ".dat").c_str()));
            ibFiles_.back() << "Title = \"" << ibObj.name() << "\"\n";

            ibForces_.push_back(std::ofstream((ibObj.name() + "Forces.dat").c_str()));
            ibForces_.back() << "Title = \"" << ibObj.name() << " Forces\"\n"
                             << "Variables = \"time\", \"x\", \"y\", \"acc_x\", \"acc_y\", \"vel_x\", \"vel_y\", \"f_norm_x\", \"f_norm_y\", \"f_shear_x\", \"f_shear_y\", \"num_fresh_cells\", \"num_dead_cells\"\n";
        }

        cutCellFile_.open("CutCells.dat");
        cutCellFile_ << "Title = \"Cut Cells\"\n";
    }
}

IbViewer::~IbViewer()
{
    close();
}

void IbViewer::close()
{
    for (std::ofstream &f: ibFiles_)
        f.close();

    for(std::ofstream &f: ibForces_)
        f.close();

    cutCellFile_.close();
}

void IbViewer::write(Scalar solutionTime, const Communicator &comm)
{
    if (solver_.comm().isMainProc())
    {
        auto geomFileItr = ibFiles_.begin();
        auto forceFileItr = ibForces_.begin();

        for (const ImmersedBoundaryObject &ibObj: ib_.ibObjs())
        {
            const Shape2D *shape = &ibObj.shape();

            //- Add a geometry record for this zone
            switch (ibObj.shape().type())
            {
            case Shape2D::CIRCLE:
            {
                const Circle &circ = *static_cast<const Circle*>(shape);

                *geomFileItr << "Geometry x=" << circ.centroid().x << " y=" << circ.centroid().y
                             << " T=CIRCLE C=BLACK LT=0.4 FC=CUST2 CS=GRID EP=300 ZN=" << zoneNo_ << "\n"
                             << circ.radius() << std::endl;
            }
                break;

            case Shape2D::BOX:
            {
                const Box& box = *static_cast<const Box*>(shape);

                *geomFileItr << "Geometry x=" << box.lower().x << " y=" << box.lower().y
                             << " T=RECTANGLE C=BLACK LT=0.3 CS=GRID ZN=" << zoneNo_ << "\n"
                             << box.upper().x - box.lower().x << " " << box.upper().y - box.lower().y << std::endl;
            }
                break;
            default:
            {
                auto pgn = shape->polygonize();

                *geomFileItr << "Geometry x=0 y=0 T=LINE C=BLACK LT=0.4 FC=CUST2 CS=GRID ZN=" << zoneNo_ << "\n"
                             << "1\n"
                             << pgn.vertices().size() << "\n";

                for (const Point2D& pt: pgn.vertices())
                    *geomFileItr << pt.x << " " << pt.y << "\n";

                geomFileItr->flush();
            }

                break;
            }

            *forceFileItr << solutionTime << " " << ibObj.position().x << " " << ibObj.position().y << " "
                          << ibObj.acceleration().x << " " << ibObj.acceleration().y << " "
                          << ibObj.velocity().x << " " << ibObj.velocity().y << " "
                          << ibObj.normalForce().x << " " << ibObj.normalForce().y << " "
                          << ibObj.shearForce().x << " " << ibObj.shearForce().y << " "
                          << ibObj.freshCells().size() << " "
                          << ibObj.deadCells().size() << std::endl;

            *geomFileItr++;
            *forceFileItr++;
        }

        for(const CutCell& cell: ib_.constructCutCells(solver_.grid().localActiveCells()))
        {
            if(cell.solid().vertices().size() == 0)
                continue;
            if(cell.isSmall())
            {
                cutCellFile_ << "Geometry x=0 y=0 T=LINE C=BLACK LT=0.2 FC=CUST2 CS=GRID ZN=" << zoneNo_ << "\n"
                             << "1\n";
            }
            else
            {
                cutCellFile_ << "Geometry x=0 y=0 T=LINE C=BLACK LT=0.2 FC=CUST1 CS=GRID ZN=" << zoneNo_ << "\n"
                             << "1\n";
            }
            cutCellFile_ << cell.solid().vertices().size() << "\n";

            for (const Point2D& pt: cell.solid().vertices())
                cutCellFile_ << pt.x << " " << pt.y << "\n";
        }

        cutCellFile_.flush();
    }

    zoneNo_++;
}
