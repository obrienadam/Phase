#include <fstream>

#include "FiniteVolume/Motion/TranslatingMotion.h"
#include "FiniteVolume/Motion/OscillatingMotion.h"
#include "FiniteVolume/Motion/SolidBodyMotion.h"
#include "FiniteVolume/Motion/MotionProfile.h"

#include "ImmersedBoundary.h"

ImmersedBoundary::ImmersedBoundary(const Input &input,
                                   const std::shared_ptr<const FiniteVolumeGrid2D> &grid,
                                   const std::shared_ptr<CellGroup> &domainCells)
    :
      grid_(grid),
      domainCells_(domainCells),
      cellStatus_(std::make_shared<FiniteVolumeField<int>>(grid, "cellStatus", FLUID_CELLS, false, false))
{
    auto ibInput = input.boundaryInput().get_child_optional("ImmersedBoundaries");

    if(ibInput)
        for (const auto &ibObjectInput: ibInput.get())
        {
            grid_->comm().printf("Initializing immersed boundary object \"%s\".\n", ibObjectInput.first.c_str());
            auto ibObj = std::make_shared<ImmersedBoundaryObject>(ibObjectInput.first);

            //- Initialize the geometry
            std::string shape = ibObjectInput.second.get<std::string>("geometry.type");
            Point2D center = Point2D(ibObjectInput.second.get<std::string>("geometry.center"));

            if (shape == "circle")
                ibObj->initCircle(center, ibObjectInput.second.get<Scalar>("geometry.radius"));
            else if (shape == "box")
                ibObj->initBox(center, ibObjectInput.second.get<Scalar>("geometry.width"),
                               ibObjectInput.second.get<Scalar>("geometry.height"));
            else if (shape == "polygon")
            {
                std::ifstream fin;
                std::vector<Point2D> verts;
                std::string filename = "case/";
                filename += ibObjectInput.second.get<std::string>("geometry.file").c_str();

                fin.open(filename.c_str());

                if (!fin.is_open())
                    throw Exception("ImmersedBoundary", "ImmersedBoundary",
                                    "failed to open file \"" + filename + "\".");

                grid_->comm().printf("Reading data for \"%s\" from file \"%s\".\n", ibObjectInput.first.c_str(),
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

                ibObj->initPolygon(verts.begin(), verts.end());
            }
            else
                throw Exception("ImmersedBoundary", "ImmersedBoundary",
                                "invalid geometry type \"" + shape + "\".");

            //- Optional geometry parameters
            boost::optional<Scalar> scaleFactor = ibObjectInput.second.get_optional<Scalar>("geometry.scale");

            if (scaleFactor)
            {
                grid_->comm().printf("Scaling \"%s\" by a factor of %lf.\n", ibObjectInput.first.c_str(),
                                     scaleFactor.get());
                ibObj->shape().scale(scaleFactor.get());
            }

            boost::optional<Scalar> rotationAngle = ibObjectInput.second.get_optional<Scalar>("geometry.rotate");

            if (rotationAngle)
            {
                grid_->comm().printf("Rotating \"%s\" by an angle of %lf degrees.\n",
                                     ibObjectInput.first.c_str(),
                                     rotationAngle.get());

                if (ibObj->shape().type() == Shape2D::BOX)
                {
                    Box *box = (Box *) &ibObj->shape();
                    auto verts = box->vertices();

                    ibObj->initPolygon(verts.begin(), verts.end());
                }

                ibObj->shape().rotate(rotationAngle.get() * M_PI / 180.);
            }

            //- Properties
            ibObj->rho = ibObjectInput.second.get<Scalar>("properties.rho", 0.);

            //- Boundary information
            for (const auto &child: ibObjectInput.second)
            {
                if (child.first == "geometry" || child.first == "interpolation"
                        || child.first == "motion" || child.first == "properties")
                    continue;

                ibObj->addBoundaryCondition(child.first, child.second.get<std::string>("type"), child.second.get<std::string>("value"));
            }

            //- Set the motion type
            std::shared_ptr<Motion> motion = nullptr;
            std::string motionType = ibObjectInput.second.get<std::string>("motion.type", "none");

            if (motionType == "translating")
            {
                motion = std::make_shared<TranslatingMotion>(
                            ibObj->position(),
                            ibObjectInput.second.get<std::string>("motion.velocity"),
                            ibObjectInput.second.get<std::string>("motion.acceleration", "(0,0)")
                            );
            }
            else if (motionType == "oscillating")
            {
                motion = std::make_shared<OscillatingMotion>(
                            ibObj->position(),
                            ibObjectInput.second.get<std::string>("motion.frequency"),
                            ibObjectInput.second.get<std::string>("motion.amplitude"),
                            ibObjectInput.second.get<std::string>("motion.phase", "(0,0)"),
                            0
                            );
            }
            else if (motionType == "solidBody" || motionType == "solid-body")
            {
                motion = std::make_shared<SolidBodyMotion>(
                            ibObj,
                            ibObjectInput.second.get<std::string>("motion.velocity", "(0,0)"),
                            (bool)ibObjectInput.second.get_child_optional("motion.axis"),
                            ibObjectInput.second.get<std::string>("motion.axis", "(0,1)")
                            );
            }
            else if(motionType == "motionProfile" || motionType == "motion-profile")
            {
                motion = std::make_shared<MotionProfile>(ibObj->position());
            }
            else if (motionType == "none")
            {
                motion = nullptr;
            }
            else
                throw Exception("ImmersedBoundary", "ImmersedBoundary", "invalid motion type \"" + motionType + "\".");

            ibObj->setMotion(motion);
            ibObjs_.push_back(ibObj);
        }

    ibInput = input.boundaryInput().get_child_optional("ImmersedBoundaryGeometryFile");

    if(ibInput)
    {
        std::string filename = ibInput.get().get<std::string>("filename");

        for(const auto &ibObjInput: input.read(filename))
        {
            std::string name = ibObjInput.first;
            auto ibObj = std::make_shared<ImmersedBoundaryObject>(name);

            ibObj->initCircle(
                        ibObjInput.second.get<std::string>("geometry.center"),
                        ibObjInput.second.get<Scalar>("geometry.radius")
                        );

            //- Init bcs, set the defaults first
            for(const auto &ibFieldBcInput: ibInput.get().get_child("fields"))
            {
                ibObj->addBoundaryCondition(ibFieldBcInput.first,
                                            ibFieldBcInput.second.get<std::string>("type"),
                                            ibFieldBcInput.second.get<std::string>("value"));
            }

            ibObjs_.push_back(ibObj);
        }
    }


    //- Collision model
    collisionModel_ = std::make_shared<CollisionModel>(
                input.boundaryInput().get<Scalar>("ImmersedBoundaries.Collisions.stiffness", 1e-4),
                input.boundaryInput().get<Scalar>("ImmersedBoundaries.Collisions.range", 0.)
                );
}

void ImmersedBoundary::setDomainCells(const std::shared_ptr<CellGroup> &domainCells)
{
    domainCells_ = domainCells;
    setCellStatus();
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

std::shared_ptr<ImmersedBoundaryObject> ImmersedBoundary::ibObj(const Point2D &pt)
{
    for (const auto &ibObj: ibObjs_)
        if (ibObj->isInIb(pt))
            return ibObj;

    return nullptr;
}

std::shared_ptr<const ImmersedBoundaryObject> ImmersedBoundary::ibObj(const Point2D &pt) const
{
    for (const auto &ibObj: ibObjs_)
        if (ibObj->isInIb(pt))
            return ibObj;

    return nullptr;
}

std::shared_ptr<ImmersedBoundaryObject> ImmersedBoundary::ibObj(const Cell &cell)
{
    return ibObj(cell.centroid());
}

std::shared_ptr<const ImmersedBoundaryObject> ImmersedBoundary::ibObj(const Cell &cell) const
{
    return ibObj(cell.centroid());
}

const std::vector<std::shared_ptr<const ImmersedBoundaryObject> > &ImmersedBoundary::findAllIbObjs(const Point2D &pt) const
{
    query_.clear();

    std::copy_if(ibObjs_.begin(), ibObjs_.end(), std::back_inserter(query_),
                 [&pt](const std::shared_ptr<ImmersedBoundaryObject> &ibObj) { return ibObj->isInIb(pt); });

    return query_;
}

std::shared_ptr<const ImmersedBoundaryObject> ImmersedBoundary::nearestIbObj(const Point2D &pt) const
{
    return nearestIntersect(pt).first;
}

std::pair<std::shared_ptr<const ImmersedBoundaryObject>, Point2D>
ImmersedBoundary::nearestIntersect(const Point2D &pt) const
{
    std::shared_ptr<const ImmersedBoundaryObject> nearestIbObj = nullptr;
    Point2D minXc;
    Scalar minDistSqr = std::numeric_limits<Scalar>::infinity();

    for (const auto &ibObj: *this)
    {
        Point2D xc = ibObj->nearestIntersect(pt);
        Scalar distSqr = (xc - pt).magSqr();

        if (distSqr < minDistSqr)
        {
            nearestIbObj = ibObj;
            minXc = xc;
            minDistSqr = distSqr;
        }
    }

    return std::make_pair(nearestIbObj, minXc);
}

std::shared_ptr<const ImmersedBoundaryObject> ImmersedBoundary::ibObj(const std::string &name) const
{
    for (const auto &ibObj: ibObjs_)
        if (ibObj->name() == name)
            return ibObj;

    throw Exception("ImmersedBoundary", "ibObj", "no immersed boundary object named \"" + name + "\".");
}

void ImmersedBoundary::updateIbPositions(Scalar timeStep)
{
    for(const auto& ibObj: ibObjs_)
        ibObj->updatePosition(timeStep);
}

FiniteVolumeEquation<Vector2D> ImmersedBoundary::velocityBcs(VectorFiniteVolumeField &u) const
{
    FiniteVolumeEquation<Vector2D> eqn(u);

    return eqn;
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

void ImmersedBoundary::applyHydrodynamicForce(Scalar rho,
                                              Scalar mu,
                                              const VectorFiniteVolumeField &u,
                                              const ScalarFiniteVolumeField &p,
                                              const Vector2D &g)
{   
    if (collisionModel_)
        for (auto ibObjP: ibObjs_)
        {
            for (auto ibObjQ: ibObjs_)
                ibObjP->applyForce(collisionModel_->force(*ibObjP, *ibObjQ));

            ibObjP->applyForce(collisionModel_->force(*ibObjP, *grid_));
        }
}

void ImmersedBoundary::applyHydrodynamicForce(const ScalarFiniteVolumeField &rho,
                                              const ScalarFiniteVolumeField &mu,
                                              const VectorFiniteVolumeField &u,
                                              const ScalarFiniteVolumeField &p,
                                              const Vector2D &g)
{   
    if (collisionModel_)
        for (auto ibObjP: ibObjs_)
        {
            for (auto ibObjQ: ibObjs_)
                ibObjP->applyForce(collisionModel_->force(*ibObjP, *ibObjQ));

            ibObjP->applyForce(collisionModel_->force(*ibObjP, *grid_));
        }
}

void ImmersedBoundary::applyCollisionForce(bool add)
{
    if(collisionModel_)
        for (auto ibObjP: ibObjs_)
        {
            //- Dont compute if no motion
            if(!ibObjP->motion())
                continue;

            Vector2D fc(0., 0.);

            //- Collisions with particles
            for (auto ibObjQ: ibObjs_)
                fc += collisionModel_->force(*ibObjP, *ibObjQ);

            if(add)
                ibObjP->addForce(fc);
            else
                ibObjP->applyForce(fc);

            //- Collisions with domain boundaries
            //            if(add)
            //                ibObjP->addForce(collisionModel_->force(*ibObjP, *grid_));
            //            else
            //                ibObjP->applyForce(collisionModel_->force(*ibObjP, *grid_));
        }
}


//- Protected

void ImmersedBoundary::setCellStatus()
{
    cellStatus_->fill(FLUID_CELLS, *domainCells_);

    for (const auto &ibObj: ibObjs_)
    {
        cellStatus_->fill(IB_CELLS, ibObj->ibCells());
        cellStatus_->fill(SOLID_CELLS, ibObj->solidCells());
    }

    grid_->sendMessages(*cellStatus_);
}
