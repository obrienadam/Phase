#include <fstream>

#include "ImmersedBoundary.h"
#include "Solver.h"
#include "GhostCellImmersedBoundary.h"

ImmersedBoundary::ImmersedBoundary(const Input &input, const Communicator& comm, Solver &solver)
    :
      solver_(solver),
      cellStatus_(solver.addScalarField("cellStatus"))
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

    for(const auto& ibObject: input.boundaryInput().get_child("ImmersedBoundaries"))
    {
        comm.printf("Initializing immersed boundary object \"%s\".\n", ibObject.first.c_str());

        //- Initialize the geometry
        const std::string type = ibObject.second.get<std::string>("geometry.type");
        const Point2D center = Point2D(ibObject.second.get<std::string>("geometry.center"));

        if(type == "circle")
        {
            ibObjs_.push_back(
                        ImmersedBoundaryObject(ibObject.first,
                                               center,
                                               ibObject.second.get<Scalar>("geometry.radius"),
                                               id++,
                                               solver.grid())
                        );
        }
        else if(type == "box")
        {
            Scalar halfWidth = ibObject.second.get<Scalar>("geometry.width")/2.,
                    halfHeight = ibObject.second.get<Scalar>("geometry.height")/2.;

            std::vector<Point2D> box = {
                center + Vector2D(-halfWidth, -halfHeight),
                center + Vector2D(halfWidth, -halfHeight),
                center + Vector2D(halfWidth, halfHeight),
                center + Vector2D(-halfWidth, halfHeight)
            };

            ibObjs_.push_back(
                        ImmersedBoundaryObject(ibObject.first,
                                               box,
                                               id++,
                                               solver.grid())
                        );
        }
        else if(type == "polygon")
        {
            std::ifstream fin;
            std::vector<Point2D> verts;
            std::string filename = "case/";
            filename += ibObject.second.get<std::string>("geometry.file").c_str();

            fin.open(filename.c_str());

            if(!fin.is_open())
                throw Exception("ImmersedBoundary", "ImmersedBoundary", "failed to open file \"" + filename + "\".");

            comm.printf("Reading data for \"%s\" from file \"%s\".\n", ibObject.first.c_str(), filename.c_str());

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

            ibObjs_.push_back(
                        ImmersedBoundaryObject(ibObject.first,
                                               verts,
                                               id++,
                                               solver.grid())
                        );
        }
        else
            throw Exception("ImmersedBoundaryObject", "ImmersedBoundaryObject", "invalid geometry type \"" + type + "\".");

        //- Optional geometry parameters
        boost::optional<Scalar> scaleFactor = ibObject.second.get_optional<Scalar>("geometry.scale");
        boost::optional<Scalar> rotationAngle = ibObject.second.get_optional<Scalar>("geometry.rotate");

        if(scaleFactor)
        {
            comm.printf("Scaling \"%s\" by a factor of %lf.\n", ibObject.first.c_str(), scaleFactor.get());
            ibObjs_.back().shape().scale(scaleFactor.get());
        }

        if(rotationAngle)
        {
            comm.printf("Rotating \"%s\" by an angle of %lf degrees.\n", ibObject.first.c_str(), rotationAngle.get());
            ibObjs_.back().shape().rotate(rotationAngle.get()*M_PI/180.);
        }

        //- Boundary information
        for(const auto& child: ibObject.second)
        {
            if(child.first == "geometry")
                continue;

            std::string type = child.second.get<std::string>("type");
            ImmersedBoundaryObject::BoundaryType boundaryType;
            Scalar boundaryRefValue = 0.;

            if(type == "fixed")
                boundaryType = ImmersedBoundaryObject::FIXED;
            else if(type == "normal_gradient")
                boundaryType = ImmersedBoundaryObject::NORMAL_GRADIENT;
            else if(type == "contact_angle")
                boundaryType = ImmersedBoundaryObject::CONTACT_ANGLE;
            else if(type == "partial_slip")
            {
                boundaryType = ImmersedBoundaryObject::PARTIAL_SLIP;
                boundaryRefValue = child.second.get<Scalar>("lambda");
            }
            else
                throw Exception("ImmersedBoundary", "ImmersedBoundary", "unrecognized boundary type \"" + type + "\".");

            comm.printf("Setting boundary type \"%s\" for field \"%s\".\n", type.c_str(), child.first.c_str());

            ibObjs_.back().addBoundaryType(child.first, boundaryType);
            ibObjs_.back().addBoundaryRefValue(child.first, boundaryRefValue);

            if(child.first == "p")
            {
                ibObjs_.back().addBoundaryType("pCorr", boundaryType);
                ibObjs_.back().addBoundaryType("dp", boundaryType);
            }
        }
    }
}

void ImmersedBoundary::initCellZones(const Communicator &comm)
{
    for(ImmersedBoundaryObject& ibObj: ibObjs_)
        ibObj.setInternalCells(comm);

    setCellStatus(comm);
}

Equation<Scalar> ImmersedBoundary::eqns(ScalarFiniteVolumeField &field)
{
    return gc::ib(ibObjs_, field);
}

Equation<Vector2D> ImmersedBoundary::eqns(VectorFiniteVolumeField &field)
{
    return gc::ib(ibObjs_, field);
}

bool ImmersedBoundary::isIbCell(const Cell &cell) const
{
    for(const ImmersedBoundaryObject& ibObj: ibObjs_)
    {
        if(ibObj.ibCells().isInGroup(cell))
            return true;
    }

    return false;
}

//- Protected

void ImmersedBoundary::setCellStatus(const Communicator &comm)
{

    for(const Cell &cell: solver_.grid().cellZone("fluid"))
        cellStatus_(cell) = FLUID;

    for(const ImmersedBoundaryObject& ibObj: ibObjs_)
    {
        for(const Cell& cell: ibObj.ibCells())
            cellStatus_(cell) = IB;

        for(const Cell& cell: ibObj.solidCells())
            cellStatus_(cell) = SOLID;
    }

    cellStatus_.grid.sendMessages(comm, cellStatus_);
}
