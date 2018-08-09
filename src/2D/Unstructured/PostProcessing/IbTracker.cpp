#include "IbTracker.h"

IbTracker::IbTracker(int fileWriteFreq,
                     const std::weak_ptr<const ImmersedBoundary> &ib,
                     double lineThickness,
                     const std::string &fillColor)
        :
        Object(fileWriteFreq),
        lineThickness_(lineThickness),
        fillColor_(fillColor),
        ib_(ib)
{
    path_ /= "IbTracker";

    if (ib_.lock()->grid()->comm().isMainProc())
    {
        createOutputDirectory();

        for (const auto &ibObj: *ib_.lock())
        {


            std::ofstream fout((path_ / (ibObj->name() + ".dat")).string());
            fout << "Title = \"" << ibObj->name() << "\"\n";
            fout.close();

            fout.open((path_ / (ibObj->name() + "_time_series.csv")).string());
            fout << "time,x,y,theta,vx,vy,omega,fx,fy,tau\n";
            fout.close();
        }
    }
}

void IbTracker::compute(Scalar time, bool force)
{
    //- Remove invalid objects
    if (do_update() || force)
    {
        //- Loop over remaining
        if (ib_.lock()->grid()->comm().isMainProc())
            for (const auto& ibObj: *ib_.lock())
            {
                std::ofstream fout((path_ / (ibObj->name() + ".dat")).string(),
                                   std::ofstream::out | std::ofstream::app);

                //- Add a geometry record for this zone
                switch (ibObj->shape().type())
                {
                    case Shape2D::CIRCLE:
                    {
                        const Circle &circ = static_cast<const Circle &>(ibObj->shape());

                        fout << "Geometry x=" << circ.centroid().x << " y=" << circ.centroid().y
                             << " T=CIRCLE C=BLACK LT=" << lineThickness_
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

                fout.open((path_ / (ibObj->name() + "_time_series.csv")).string(),
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
