#ifndef PHASE_EQUATION_H
#define PHASE_EQUATION_H

#include "SparseMatrixSolver.h"

class Equation
{

public:

    Equation(Size rank, Size nnz);

    Equation(const Equation &eqn) = default;

    Equation(Equation &&eqn) = default;

    Equation &operator=(const Equation &eqn);

    Equation &operator=(Equation &&eqn);

    Scalar coeff(Index localRow, Index globalCol) const;

    Scalar &coeffRef(Index localRow, Index globalCol);

    void clear();

    Equation &operator+=(const Equation &rhs);

    Equation &operator-=(const Equation &rhs);

    Equation &operator+=(const Vector &rhs);

    Equation &operator-=(const Vector &rhs);

    Equation &operator*=(Scalar rhs);

    Equation &operator==(Scalar rhs);

    Equation &operator==(const Equation &rhs);

    Equation &operator==(const Vector &rhs);

protected:

    void addCoeff(Index localRow, Index globalCol, Scalar val);

    void setCoeff(Index localRow, Index globalCol, Scalar val);

    void addRhs(Index localRow, Scalar val)
    { _rhs(localRow) += val; }

    void setRhs(Index localRow, Scalar val)
    { _rhs(localRow) = val; }

    SparseMatrixSolver::CoefficientList _coeffs;

    Vector _rhs;

    std::shared_ptr<SparseMatrixSolver> _spSolver;

};

Equation operator+(Equation lhs, const Equation &rhs);

Equation operator+(Equation lhs, const Vector &rhs);

Equation operator-(Equation lhs, const Equation &rhs);

Equation operator-(Equation lhs, const Vector &rhs);

Equation operator*(Equation lhs, Scalar rhs);

Equation operator*(Scalar lhs, Equation rhs);

#endif
