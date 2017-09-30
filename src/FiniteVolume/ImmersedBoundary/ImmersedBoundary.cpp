#include <fstream>

#include "ImmersedBoundary.h"
#include "SurfaceTensionForce.h"
#include "StepImmersedBoundaryObject.h"
#include "GhostCellImmersedBoundaryObject.h"
#include "TranslatingMotion.h"
#include "OscillatingMotion.h"
#include "SolidBodyMotion.h"

ImmersedBoundary::ImmersedBoundary(const Input &input, Solver &solver)
        :
        solver_(solver),
        cellStatus_(solver.addIntegerField("cellStatus"))
{
    if (input.boundaryInput().find("ImmersedBoundaries") == input.boundaryInput().not_found())
    {
        solver.grid().comm().printf("No immersed boundaries present.\n");
        return;
    }

    Label id = 0;

    for (const auto &ibObjectInput: input.boundaryInput().get_child("ImmersedBoundaries"))
    {
        solver.grid().comm().printf("Initializing immersed boundary object \"%s\".\n", ibObjectInput.first.c_str());

        std::string method = ibObjectInput.second.get<std::string>("method", "ghost-cell");
        std::transform(method.begin(), method.end(), method.begin(), ::tolower);
        solver.grid().comm().printf("Immersed boundary method: %s\n", method.c_str());

        std::shared_ptr<ImmersedBoundaryObject> ibObject;

        if (method == "step")
            ibObject = std::make_shared<StepImmersedBoundaryObject>(ibObjectInput.first, id++, solver.grid());
        else if (method == "ghost-cell")
            ibObject = std::make_shared<GhostCellImmersedBoundaryObject>(ibObjectInput.first, id++, solver.grid());
        else
            throw Exception("ImmersedBoundary", "ImmersedBoundary",
                            "invalid immersed boundary method \"" + method + "\".");

        //- Initialize the geometry
        std::string shape = ibObjectInput.second.get<std::string>("geometry.type");
        Point2D center = Point2D(ibObjectInput.second.get<std::string>("geometry.center"));

        if (shape == "circle")
        {
            ibObject->initCircle(center, ibObjectInput.second.get<Scalar>("geometry.radius"));
        }
        else if (shape == "box")
        {
            ibObject->initBox(center, ibObjectInput.second.get<Scalar>("geometry.width"),
                              ibObjectInput.second.get<Scalar>("geometry.height"));
        }
        else if (shape == "polygon")
        {
            std::ifstream fin;
            std::vector<Point2D> verts;
            std::string filename = "case/";
            filename += ibObjectInput.second.get<std::string>("geometry.file").c_str();

            fin.open(filename.c_str());

            if (!fin.is_open())
                throw Exception("ImmersedBoundary", "ImmersedBoundary", "failed to open file \"" + filename + "\".");

            solver.grid().comm().printf("Reading data for \"%s\" from file \"%s\".\n", ibObjectInput.first.c_str(),
                                        filename.c_str());

            while (!fin.eof())
            {
                Scalar x, y;
                fin >> x;
                fin >> y;

                verts.push_back(Point2D(x, y));
            }

            fin.close();

            Vector2D translation = center - Polygon(verts.begin(), verts.end()).centroid();

            for (Point2D &vert: verts)
                vert += translation;

            ibObject->initPolygon(verts.begin(), verts.end());
        }
        else
            throw Exception("ImmersedBoundaryObject", "ImmersedBoundaryObject",
                            "invalid geometry type \"" + shape + "\".");

        //- Optional geometry parameters
        boost::optional<Scalar> scaleFactor = ibObjectInput.second.get_optional<Scalar>("geometry.scale");

        if (scaleFactor)
        {
            solver.grid().comm().printf("Scaling \"%s\" by a factor of %lf.\n", ibObjectInput.first.c_str(),
                                        scaleFactor.get());
            ibObject->shape().scale(scaleFactor.get());
        }

        boost::optional<Scalar> rotationAngle = ibObjectInput.second.get_optional<Scalar>("geometry.rotate");

        if (rotationAngle)
        {
            solver.grid().comm().printf("Rotating \"%s\" by an angle of %lf degrees.\n", ibObjectInput.first.c_str(),
                                        rotationAngle.get());

            if (ibObject->shape().type() == Shape2D::BOX)
            {
                Box *box = (Box *) &ibObject->shape();
                auto verts = box->vertices();

                ibObject->initPolygon(verts.begin(), verts.end());
            }

            ibObject->shape().rotate(rotationAngle.get() * M_PI / 180.);
        }

        //- Properties
        ibObject->rho = ibObjectInput.second.get<Scalar>("properties.rho", 0.);

        //- Boundary information
        for (const auto &child: ibObjectInput.second)
        {
            if (child.first == "geometry" || child.first == "interpolation"
                || child.first == "motion" || child.first == "method" || child.first == "properties")
                continue;

            std::string type = child.second.get<std::string>("type");
            ImmersedBoundaryObject::BoundaryType boundaryType;

            if (type == "fixed")
            {
                boundaryType = ImmersedBoundaryObject::FIXED;
                ibObject->addBoundaryRefValue(child.first, child.second.get<std::string>("value"));
            }
            else if (type == "normal_gradient")
            {
                boundaryType = ImmersedBoundaryObject::NORMAL_GRADIENT;
                ibObject->addBoundaryRefValue(child.first, child.second.get<std::string>("value"));
            }
            else if (type == "partial_slip")
            {
                boundaryType = ImmersedBoundaryObject::PARTIAL_SLIP;
                ibObject->addBoundaryRefValue(child.first, child.second.get<std::string>("value"));
            }
            else
                throw Exception("ImmersedBoundary", "ImmersedBoundary", "unrecognized boundary type \"" + type + "\".");

            solver.grid().comm().printf("Setting boundary type \"%s\" for field \"%s\".\n", type.c_str(),
                                        child.first.c_str());
            ibObject->addBoundaryType(child.first, boundaryType);
        }

        //- Set the motion type
        std::shared_ptr<Motion> motion = nullptr;
        std::string motionType = ibObjectInput.second.get<std::string>("motion.type", "none");

        if (motionType == "translating")
        {
            motion = std::make_shared<TranslatingMotion>(
                    ibObject,
                    ibObjectInput.second.get<std::string>("motion.velocity"),
                    ibObjectInput.second.get<std::string>("motion.acceleration", "(0,0)")
            );
        }
        else if (motionType == "oscillating")
        {
            motion = std::make_shared<OscillatingMotion>(
                    ibObject,
                    ibObjectInput.second.get<std::string>("motion.frequency"),
                    ibObjectInput.second.get<std::string>("motion.amplitude"),
                    ibObjectInput.second.get<std::string>("motion.phase", "(0,0)"),
                    0
            );
        }
        else if (motionType == "solidBody")
        {
            motion = std::make_shared<SolidBodyMotion>(ibObject);
        }
        else if (motionType == "none")
        {
            motion = nullptr;
        }
        else
            throw Exception("ImmersedBoundary", "ImmersedBoundary", "invalid motion type \"" + motionType + "\".");

        ibObject->setMotion(motion);

        ibObjs_.push_back(ibObject);
    }
}

void ImmersedBoundary::initCellZones(CellZone &zone)
{
    zone_ = &zone;

    for (auto &ibObj: ibObjs_)
    {
        ibObj->setZone(zone);
        ibObj->updateCells();
    }

    setCellStatus();
    solver_.grid().computeGlobalOrdering();
}

void ImmersedBoundary::clearCellZones()
{
    for (auto &ibObj: ibObjs_)
    {
        ibObj->clear();
    }

    solver_.grid().computeGlobalOrdering();
}

CellGroup ImmersedBoundary::ibCells() const
{
    CellGroup ibCellGroup;

    for (const auto &ibObj: ibObjs_)
        ibCellGroup += ibObj->ibCells();

    return ibCellGroup;
}

CellGroup ImmersedBoundary::solidCells() const
{
    CellGroup solidCellGroup;

    for (const auto &ibObj: ibObjs_)
        solidCellGroup += ibObj->solidCells();

    return solidCellGroup;
}

std::shared_ptr<const ImmersedBoundaryObject> ImmersedBoundary::ibObj(const Point2D &pt) const
{
    for (const auto &ibObj: ibObjs_)
        if (ibObj->isInIb(pt))
            return ibObj;

    return nullptr;
}

const ImmersedBoundaryObject &ImmersedBoundary::ibObj(const std::string &name) const
{
    for (const auto &ibObj: ibObjs_)
        if (ibObj->name() == name)
            return *ibObj;

    throw Exception("ImmersedBoundary", "ibObj", "no immersed boundary object named \"" + name + "\".");
}

std::vector<Ref<const ImmersedBoundaryObject> > ImmersedBoundary::ibObjs() const
{
    std::vector<Ref<const ImmersedBoundaryObject>> refs;

    std::transform(ibObjs_.begin(), ibObjs_.end(), std::back_inserter(refs),
                   [](const std::shared_ptr<ImmersedBoundaryObject> &ptr) { return std::cref(*ptr); });

    return refs;
}

void ImmersedBoundary::update(Scalar timeStep)
{
    for (auto &ibObj: ibObjs_)
        ibObj->update(timeStep);

    setCellStatus();
    solver_.grid().computeGlobalOrdering();
}

Equation<Vector2D> ImmersedBoundary::velocityBcs(VectorFiniteVolumeField &u) const
{
    Equation<Vector2D> eqn(u);

    for (const auto &ibObj: ibObjs_)
        eqn += ibObj->velocityBcs(u);

    return eqn;
}

Equation<Scalar> ImmersedBoundary::pressureBcs(Scalar rho, ScalarFiniteVolumeField &p) const
{
    Equation<Scalar> eqn(p);

    for (const auto &ibObj: ibObjs_)
        eqn += ibObj->pressureBcs(rho, p);

    return eqn;
}

Equation<Scalar> ImmersedBoundary::contactLineBcs(const SurfaceTensionForce &fst, ScalarFiniteVolumeField &gamma) const
{
    Equation<Scalar> eqn(gamma);

    for (const auto &ibObj: ibObjs_)
        eqn += ibObj->contactLineBcs(gamma, fst.getTheta(*ibObj));

    return eqn;
}

void ImmersedBoundary::clearFreshCells()
{
    for (auto &ibObj: ibObjs_)
        ibObj->clearFreshCells();
}

bool ImmersedBoundary::isIbCell(const Cell &cell) const
{
    for (const auto &ibObj: ibObjs_)
    {
        if (ibObj->isInIb(cell.centroid()))
            return true;
    }

    return false;
}

void ImmersedBoundary::computeForce(Scalar rho,
                                    Scalar mu,
                                    const VectorFiniteVolumeField &u,
                                    const ScalarFiniteVolumeField &p,
                                    const Vector2D &g)
{
    for (auto ibObj: ibObjs_)
        ibObj->computeForce(rho, mu, u, p, g);

    if(collisionModel_)
        for (auto ibObjP: ibObjs_)
            for(auto ibObjQ: ibObjs_)
            {
                if(ibObjP == ibObjQ)
                    continue;

                collisionModel_->force(*ibObjP, *ibObjQ);
            }
}

void ImmersedBoundary::computeForce(const ScalarFiniteVolumeField &rho,
                                    const ScalarFiniteVolumeField &mu,
                                    const VectorFiniteVolumeField &u,
                                    const ScalarFiniteVolumeField &p)
{

}

//- Protected

void ImmersedBoundary::setCellStatus()
{
    for (const Cell &cell: solver_.grid().cellZone("fluid"))
        cellStatus_(cell) = FLUID_CELLS;

    for (const CellZone &bufferZone: solver_.grid().bufferZones())
        for (const Cell &cell: bufferZone)
            cellStatus_(cell) = BUFFER_CELLS;

    for (const auto &ibObj: ibObjs_)
    {
        for (const Cell &cell: ibObj->ibCells())
            cellStatus_(cell) = IB_CELLS;

        for (const Cell &cell: ibObj->solidCells())
            cellStatus_(cell) = SOLID_CELLS;

        for (const Cell &cell: ibObj->freshCells())
            cellStatus_(cell) = FRESH_CELLS;

        for (const Cell &cell: ibObj->deadCells())
            cellStatus_(cell) = DEAD_CELLS;
    }
}
