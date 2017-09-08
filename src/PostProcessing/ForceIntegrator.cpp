#include "ForceIntegrator.h"

ForceIntegrator::ForceIntegrator(const Solver &solver,
                                 const Patch &patch,
                                 const VectorFiniteVolumeField &u,
                                 const ScalarFiniteVolumeField &rho,
                                 const ScalarFiniteVolumeField &mu,
                                 const ScalarFiniteVolumeField &p)
        :
        PostProcessingObject(solver),
        patch_(patch),
        u_(u),
        rho_(rho),
        mu_(mu),
        p_(p)
{

}

void ForceIntegrator::compute(Scalar time)
{
    Vector2D f(0., 0.);

    for (const Face &face: patch_)
    {
        Vector2D sf = face.outwardNorm();
        Vector2D tf = sf.tangentVec();
        Scalar l = (face.lCell().centroid() - face.centroid()).mag();

        f += p_(face) * sf + 0.5 * rho_(face) * pow(dot(u_(face), sf), 2) / sf.magSqr()
             + mu_(face) * dot(u_(face.lCell()) - u_(face), tf) / l;
    }

    time_.push_back(time);
    force_.push_back(f);
}

void ForceIntegrator::write() const
{
    std::ofstream fout;
    fout.open(patch_.name() + "_force.dat");

    fout << "title = \"" << patch_.name() << " force\"\n"
         << "variables = \"time\", \"force_x\", \"force_y\"\n"
         << "zone I = " << time_.size() << ", F=BLOCK\n";

    for (Scalar time: time_)
        fout << time << "\n";

    for (const Vector2D &f: force_)
        fout << f.x << "\n";

    for (const Vector2D &f: force_)
        fout << f.y << "\n";

    fout.close();
}
