#include "IbTracker.h"
#include "Solver.h"

IbTracker::IbTracker(const Solver &solver, double lineThickness, const std::string &fillColor)
        :
        PostProcessingObject(solver),
        lineThickness_(lineThickness),
        fillColor_(fillColor)
{
    outputDir_ = outputDir_ / "IbTracker";

    if (solver.grid()->comm().isMainProc())
    {
        createOutputDirectory();

        for (auto ibObj: solver.ib())
        {
            std::cout << ibObj->name() << "\n";
            std::ofstream fout((outputDir_ / (ibObj->name() + ".dat")).string());
            fout << "Title = \"" << ibObj->name() << "\"\n";
            fout.close();

            fout.open((outputDir_ / (ibObj->name() + "_timeSeriesPlot.dat")).string());
            fout.close();
        }
    }

    ibObjs_.insert(ibObjs_.end(), solver_.ib().begin(), solver_.ib().end());
}

void IbTracker::compute(Scalar time)
{
    //- Remove invalid objects
    if (iterNo_++ % fileWriteFrequency_ == 0)
    {
        ibObjs_.erase(
                std::remove_if(
                        ibObjs_.begin(),
                        ibObjs_.end(),
                        [](std::weak_ptr<ImmersedBoundaryObject> ibObj) {
                            return !ibObj.lock();
                        }),
                ibObjs_.end());

        //- Loop over remaining
        if (solver_.grid()->comm().isMainProc())
            for (auto ptr: ibObjs_)
            {
                auto ibObj = ptr.lock();

                std::ofstream fout((outputDir_ / (ibObj->name() + ".dat")).string(),
                                   std::ofstream::out | std::ofstream::app);

                //- Add a geometry record for this zone
                switch (ibObj->shape().type())
                {
                    case Shape2D::CIRCLE:
                    {
                        const Circle &circ = static_cast<const Circle &>(ibObj->shape());

                        fout << "Geometry x=" << circ.centroid().x << " y=" << circ.centroid().y
                             << " T=CIRCLE C=BLACK LT=" << lineThickness_ << " FC=" << fillColor_
                             << " CS=GRID EP=300 ZN=" << zoneNo_ << "\n" << circ.radius() << std::endl;
                    }
                        break;

                    case Shape2D::BOX:
                    {
                        const Box &box = static_cast<const Box &>(ibObj->shape());

                        fout << "Geometry x=" << box.lower().x << " y=" << box.lower().y
                             << " T=RECTANGLE C=BLACK LT=" << lineThickness_ << " FC=" << fillColor_
                             << " CS=GRID ZN=" << zoneNo_ << "\n"
                             << box.upper().x - box.lower().x << " " << box.upper().y - box.lower().y << std::endl;
                    }
                        break;
                    default:
                    {
                        auto pgn = ibObj->shape().polygonize();

                        fout << "Geometry x=0 y=0 T=LINE C=BLACK LT=" << lineThickness_ << " FC=" << fillColor_
                             << " CS=GRID ZN=" << zoneNo_ << "\n" << "1\n" << pgn.vertices().size() << "\n";

                        for (const Point2D &pt: pgn.vertices())
                            fout << pt.x << " " << pt.y << "\n";
                    }
                }

                fout.close();

                fout.open((outputDir_ / (ibObj->name() + "_timeSeriesPlot.dat")).string(),
                          std::ofstream::out | std::ofstream::app);

                fout << time << ","
                     << ibObj->position().x << ","
                     << ibObj->position().y << ","
                     << ibObj->theta() << ","
                     << ibObj->velocity().x << ","
                     << ibObj->velocity().y << ","
                     << ibObj->omega() << ","
                     << ibObj->force().x << ","
                     << ibObj->force().y << ","
                     << ibObj->torque() << "\n";

                fout.close();
            }

        zoneNo_++;
    }
}