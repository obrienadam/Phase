#include "JacobianField.h"

JacobianField::JacobianField(const VectorFiniteVolumeField& u, const std::shared_ptr<CellGroup> &cells)
        :
        TensorFiniteVolumeField(u.grid(),
                                u.name() + "Jacobian",
                                Tensor2D(),
                                true,
                                false,
                                cells),
        u_(u)
{

}

void JacobianField::computeFaces()
{
    for(const Face& f: grid_->interiorFaces())
    {
        Vector2D rc = f.rCell().centroid() - f.lCell().centroid();
        (*this)(f) = outer(u_(f.rCell()) - u_(f.lCell()), rc/rc.magSqr());
    }

    for(const Face& f: grid_->boundaryFaces())
    {
        Vector2D rf = f.centroid() - f.lCell().centroid();
        (*this)(f) = outer(u_(f) - u_(f.lCell()), rf/rf.magSqr());
    }
}

void JacobianField::compute(const CellGroup& cells)
{
    computeFaces();

    for(const Cell& cell: cells)
    {
        Vector2D sumA(0., 0.);
        Tensor2D tmp(0., 0., 0., 0.);

        for(const InteriorLink& nb: cell.neighbours())
        {
            const Tensor2D &jf = (*this)(nb.face());
            Vector2D sf = nb.outwardNorm().abs();

            tmp += Tensor2D(
                    jf.xx*sf.x, jf.xy*sf.y,
                    jf.yx*sf.x, jf.yy*sf.y
            );

            sumA += sf;
        }

        for(const BoundaryLink& bd: cell.boundaries())
        {
            const Tensor2D &jf = (*this)(bd.face());
            Vector2D sf = bd.outwardNorm().abs();

            tmp += Tensor2D(
                    jf.xx*sf.x, jf.xy*sf.y,
                    jf.yx*sf.x, jf.yy*sf.y
            );

            sumA += sf;
        }

        (*this)(cell) = Tensor2D(
                tmp.xx/sumA.x, tmp.xy/sumA.y,
                tmp.yx/sumA.x, tmp.yy/sumA.y
        );
    }
}

void JacobianField::compute()
{
    compute(*cellGroup_);
}
