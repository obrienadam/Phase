#include <fstream>

#include "ImmersedBoundary.h"
#include "TranslatingImmersedBoundaryObject.h"
#include "OscillatingImmersedBoundaryObject.h"
#include "Solver.h"
#include "GhostCellImmersedBoundary.h"

ImmersedBoundary::ImmersedBoundary(const Input &input, const Communicator& comm, Solver &solver)
    :
      solver_(solver),
      comm_(comm),
      cellStatus_(solver.addIntegerField("cellStatus"))
{
    try //- Lazy way to check if any immersed boundary input is present, just catch the exception if it fails
    {
        input.boundaryInput().get_child("ImmersedBoundaries");
    }
    catch(...)
    {
        comm.printf("No immersed boundaries present.\n");
        return;
    }

    Label id = 0;

    for(const auto& ibObjectInput: input.boundaryInput().get_child("ImmersedBoundaries"))
    {
        comm.printf("Initializing immersed boundary object \"%s\".\n", ibObjectInput.first.c_str());

        std::string motion = ibObjectInput.second.get<std::string>("motion.type", "none");

        comm.printf("Immersed boundary motion type: %s\n", motion.c_str());

        std::shared_ptr<ImmersedBoundaryObject> ibObject;

        if(motion == "none")
            ibObject = std::make_shared<ImmersedBoundaryObject>(ImmersedBoundaryObject(ibObjectInput.first, id++, solver.grid()));
        else if(motion == "translating")
        {
            ibObject = std::make_shared<TranslatingImmersedBoundaryObject>(
                        TranslatingImmersedBoundaryObject(ibObjectInput.first,
                                                          ibObjectInput.second.get<std::string>("motion.velocity", "(0,0)"),
                                                          ibObjectInput.second.get<Scalar>("motion.omega", 0),
                                                          id++,
                                                          solver.grid())
                        );
        }
        else if(motion == "oscillating")
        {
            ibObject = std::make_shared<OscillatingImmersedBoundaryObject>(
                        OscillatingImmersedBoundaryObject(ibObjectInput.first,
                                                          ibObjectInput.second.get<Scalar>("motion.amplitude", 1.),
                                                          ibObjectInput.second.get<Scalar>("motion.frequency", 0.),
                                                          id++,
                                                          solver.grid())
                        );
        }
        else
            throw Exception("ImmersedBoundary", "ImmersedBoundary", "invalid motion type + \"" + motion + "\".");

        //- Initialize the geometry
        std::string shape = ibObjectInput.second.get<std::string>("geometry.type");
        Point2D center = Point2D(ibObjectInput.second.get<std::string>("geometry.center"));

        if(shape == "circle")
        {
            ibObject->initCircle(center, ibObjectInput.second.get<Scalar>("geometry.radius"));
        }
        else if(shape == "box")
        {
            Scalar halfWidth = ibObjectInput.second.get<Scalar>("geometry.width")/2.,
                    halfHeight = ibObjectInput.second.get<Scalar>("geometry.height")/2.;

            std::vector<Point2D> box = {
                center + Vector2D(-halfWidth, -halfHeight),
                center + Vector2D(halfWidth, -halfHeight),
                center + Vector2D(halfWidth, halfHeight),
                center + Vector2D(-halfWidth, halfHeight)
            };

            ibObject->initPolygon(box);
        }
        else if(shape == "polygon")
        {
            std::ifstream fin;
            std::vector<Point2D> verts;
            std::string filename = "case/";
            filename += ibObjectInput.second.get<std::string>("geometry.file").c_str();

            fin.open(filename.c_str());

            if(!fin.is_open())
                throw Exception("ImmersedBoundary", "ImmersedBoundary", "failed to open file \"" + filename + "\".");

            comm.printf("Reading data for \"%s\" from file \"%s\".\n", ibObjectInput.first.c_str(), filename.c_str());

            while(!fin.eof())
            {
                Scalar x, y;
                fin >> x;
                fin >> y;

                verts.push_back(Point2D(x, y));
            }

            fin.close();

            Vector2D translation = center - Polygon(verts).centroid();

            for(Point2D &vert: verts)
                vert += translation;

            ibObject->initPolygon(verts);
        }
        else
            throw Exception("ImmersedBoundaryObject", "ImmersedBoundaryObject", "invalid geometry type \"" + shape + "\".");

        //- Optional geometry parameters
        boost::optional<Scalar> scaleFactor = ibObjectInput.second.get_optional<Scalar>("geometry.scale");
        boost::optional<Scalar> rotationAngle = ibObjectInput.second.get_optional<Scalar>("geometry.rotate");

        if(scaleFactor)
        {
            comm.printf("Scaling \"%s\" by a factor of %lf.\n", ibObjectInput.first.c_str(), scaleFactor.get());
            ibObject->shape().scale(scaleFactor.get());
        }

        if(rotationAngle)
        {
            comm.printf("Rotating \"%s\" by an angle of %lf degrees.\n", ibObjectInput.first.c_str(), rotationAngle.get());
            ibObject->shape().rotate(rotationAngle.get()*M_PI/180.);
        }

        //- Boundary information
        for(const auto& child: ibObjectInput.second)
        {
            if(child.first == "geometry" || child.first == "interpolation" || child.first == "motion")
                continue;

            std::string type = child.second.get<std::string>("type");
            ImmersedBoundaryObject::BoundaryType boundaryType;

            if(type == "fixed")
            {
                boundaryType = ImmersedBoundaryObject::FIXED;
                ibObject->addBoundaryRefValue(child.first, child.second.get<std::string>("value"));
            }
            else if(type == "normal_gradient")
            {
                boundaryType = ImmersedBoundaryObject::NORMAL_GRADIENT;
                ibObject->addBoundaryRefValue(child.first, child.second.get<std::string>("value"));
            }
            else if(type == "contact_angle")
            {
                boundaryType = ImmersedBoundaryObject::CONTACT_ANGLE;
                ibObject->addBoundaryRefValue(child.first, child.second.get<std::string>("value"));
            }
            else if(type == "partial_slip")
            {
                boundaryType = ImmersedBoundaryObject::PARTIAL_SLIP;
                ibObject->addBoundaryRefValue(child.first, child.second.get<Scalar>("lambda"));
            }
            else
                throw Exception("ImmersedBoundary", "ImmersedBoundary", "unrecognized boundary type \"" + type + "\".");

            comm.printf("Setting boundary type \"%s\" for field \"%s\".\n", type.c_str(), child.first.c_str());

            ibObject->addBoundaryType(child.first, boundaryType);

            if(child.first == "p")
            {
                ibObject->addBoundaryType("pCorr", boundaryType);
                ibObject->addBoundaryType("dp", boundaryType);
            }
        }

        ibObjs_.push_back(ibObject);
    }
}

void ImmersedBoundary::initCellZones()
{
    for(auto& ibObj: ibObjs_)
        ibObj->setInternalCells();

    setCellStatus();
    solver_.grid().computeGlobalOrdering(solver_.comm());
}

void ImmersedBoundary::update(Scalar timeStep)
{
    for(auto& ibObj: ibObjs_)
        ibObj->update(timeStep);

    setCellStatus();
    solver_.grid().computeGlobalOrdering(comm_);
}

std::vector<Ref<const ImmersedBoundaryObject> > ImmersedBoundary::ibObjs() const
{
    std::vector<Ref<const ImmersedBoundaryObject>> refs;

    std::transform(ibObjs_.begin(), ibObjs_.end(), std::back_inserter(refs),
                   [](const std::shared_ptr<ImmersedBoundaryObject>& ptr)->const ImmersedBoundaryObject& { return *ptr; });

    return refs;
}

bool ImmersedBoundary::isIbCell(const Cell &cell) const
{
    for(const auto& ibObj: ibObjs_)
    {
        if(ibObj->ibCells().isInGroup(cell))
            return true;
    }

    return false;
}

//- Protected

void ImmersedBoundary::setCellStatus()
{

    for(const Cell &cell: solver_.grid().cellZone("fluid"))
        cellStatus_(cell) = FLUID;

    for(const auto& ibObj: ibObjs_)
    {
        for(const Cell& cell: ibObj->ibCells())
            cellStatus_(cell) = IB;

        for(const Cell& cell: ibObj->solidCells())
            cellStatus_(cell) = SOLID;

        for(const Cell& cell: ibObj->freshlyClearedCells())
            cellStatus_(cell) = FRESHLY_CLEARED;
    }

    cellStatus_.grid.sendMessages(comm_, cellStatus_);
}
