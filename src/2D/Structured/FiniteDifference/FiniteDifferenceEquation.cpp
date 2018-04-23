#include "FiniteDifferenceEquation.h"

template<>
FiniteDifferenceEquation<Scalar>::FiniteDifferenceEquation(const std::shared_ptr<FiniteDifferenceField<Scalar> > &field)
    :
      Equation(field->grid()->nNodes(), 5),
      _field(field)
{

}

template<>
FiniteDifferenceEquation<Vector2D>::FiniteDifferenceEquation(const std::shared_ptr<FiniteDifferenceField<Vector2D> > &field)
    :
      Equation(2 * field->grid()->nNodes(), 5),
      _field(field)
{

}

template<>
void FiniteDifferenceEquation<Scalar>::add(const StructuredGrid2D::Node &cnode,
                                           const StructuredGrid2D::Node &nnode,
                                           Scalar coeff)
{
    addCoeff(_field->idxMap()->local(cnode, 0), _field->idxMap()->global(nnode, 0), coeff);
}

template<>
void FiniteDifferenceEquation<Vector2D>::add(const StructuredGrid2D::Node &cnode,
                                             const StructuredGrid2D::Node &nnode,
                                             Scalar coeff)
{
    addCoeff(_field->idxMap()->local(cnode, 0), _field->idxMap()->global(nnode, 0), coeff);
    addCoeff(_field->idxMap()->local(cnode, 1), _field->idxMap()->global(nnode, 1), coeff);
}

template<>
void FiniteDifferenceEquation<Vector2D>::add(const StructuredGrid2D::Node &cnode,
                                             const StructuredGrid2D::Node &nnode,
                                             const Vector2D &coeff)
{
    addCoeff(_field->idxMap()->local(cnode, 0), _field->idxMap()->global(nnode, 0), coeff.x);
    addCoeff(_field->idxMap()->local(cnode, 1), _field->idxMap()->global(nnode, 1), coeff.y);
}
