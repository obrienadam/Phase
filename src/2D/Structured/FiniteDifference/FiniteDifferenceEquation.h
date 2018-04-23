#ifndef PHASE_FINITE_DIFFERENCE_EQUATION_H
#define PHASE_FINITE_DIFFERENCE_EQUATION_H

#include <vector>

#include "System/Input.h"
#include "Math/Equation.h"
#include "FiniteDifferenceField.h"


template<class T>
class FiniteDifferenceEquation: public Equation
{
public:

    FiniteDifferenceEquation(const std::shared_ptr<FiniteDifferenceField<T>> &field);

    FiniteDifferenceEquation(const Input &input, const std::shared_ptr<FiniteDifferenceField<T>> &field);

    void add(const StructuredGrid2D::Node &cnode, const StructuredGrid2D::Node &nnode, Scalar coeff);

    void add(const StructuredGrid2D::Node &cnode, const StructuredGrid2D::Node &nnode, const T &coeff);

    void addRhs(const StructuredGrid2D::Node &node, const T &coeff);

protected:

    std::shared_ptr<FiniteDifferenceField<T>> _field;
};

#include "FiniteDifferenceEquation.tpp"

#endif
