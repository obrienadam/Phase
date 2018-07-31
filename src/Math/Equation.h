#ifndef PHASE_EQUATION_H
#define PHASE_EQUATION_H

#include "SparseMatrixSolver.h"

class Equation
{

public:

    Equation() {}

    Equation(Size rank, Size nnz);

    Equation(const Equation &eqn) = default;

    Equation(Equation &&eqn) = default;

    Equation &operator=(const Equation &eqn);

    Equation &operator=(Equation &&eqn);

    Size rank() const
    { return _coeffs.size(); }

    void setRank(Size rank);

    void setRank(Size nRows, Size nCols);

    Scalar coeff(Index localRow, Index globalCol) const;

    Scalar &coeffRef(Index localRow, Index globalCol);

    void setCoeff(Index localRow, Index globalCol, Scalar val);

    void addCoeff(Index localRow, Index globalCol, Scalar val);

    const SparseMatrixSolver::CoefficientList &coeffs() const
    { return _coeffs; }

    const Vector &rhs() const
    { return _rhs; }

    Scalar x(Index idx) const
    { return _spSolver->x(idx); }

    Scalar b(Index idx) const
    { return _rhs(idx); }

    template<class ColIt, class ValIt>
    void setCoeffs(Index row, ColIt colBegin, ColIt colEnd, ValIt valBegin)
    {
        while(colBegin != colEnd)
            setCoeff(row, *(colBegin++), *(valBegin++));
    }

    template<class ColIt, class ValIt>
    void addCoeffs(Index row, ColIt colBegin, ColIt colEnd, ValIt valBegin)
    {
        while(colBegin != colEnd)
            addCoeff(row, *(colBegin++), *(valBegin++));
    }

    void setCoeffs(Index row, const std::initializer_list<Index> &cols, const std::initializer_list<Scalar> &vals)
    { setCoeffs(row, cols.begin(), cols.end(), vals.begin()); }

    void addCoeffs(Index row, const std::initializer_list<Index> &cols, const std::initializer_list<Scalar> &vals)
    { addCoeffs(row, cols.begin(), cols.end(), vals.begin()); }

    void addRhs(Index localRow, Scalar val)
    { _rhs(localRow) += val; }

    void setRhs(Index localRow, Scalar val)
    { _rhs(localRow) = val; }

    void setSparseSolver(const std::shared_ptr<SparseMatrixSolver> &spSolver)
    { _spSolver = spSolver; }

    virtual Scalar solve();

    virtual Scalar solveLeastSquares();

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
