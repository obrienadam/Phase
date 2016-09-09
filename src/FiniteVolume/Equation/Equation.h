#ifndef EQUATION_H
#define EQUATION_H

#include "SparseMatrix.h"
#include "FiniteVolumeGrid2D.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"

template<class T>
class Equation
{
public:

    typedef Eigen::Triplet<Scalar> Triplet;

    Equation(T& field, const std::string& name = "Unknown");
    Equation(const Input& input, T& field, const std::string& name = "Unknown");

    Equation(const Equation<T>& other);

    SparseMatrix& matrix(){ return spMat_; }
    SparseVector& boundaries(){ return boundaries_; }
    SparseVector& sources(){ return sources_; }

    Equation<T>& operator+=(const Equation<T>& rhs);
    Equation<T>& operator-=(const Equation<T>& rhs);
    Equation<T>& operator+=(const T& rhs);
    Equation<T>& operator-=(const T& rhs);

    Equation<T>& operator=(const Equation<T>& rhs);

    Equation<T>& operator*=(Scalar rhs);
    Equation<T>& operator*=(const ScalarFiniteVolumeField& rhs);

    Equation<T>& operator==(Scalar rhs);
    Equation<T>& operator==(const Equation<T>& rhs);
    Equation<T>& operator==(const T& rhs);

    Equation<T>& assemble(const std::vector<Triplet> &triplets);
    void computeSolver() const { solver_.compute(spMat_); }

    Scalar solve(bool recomputePreconditioner = true);
    Scalar solve(const SparseVector& x0, bool recomputePreconditioner = true);

    Scalar error() const { return solver_.error(); }
    int iterations() const { return solver_.iterations(); }

    const SparseMatrix& matrix() const { return spMat_; }

    void relax(Scalar relaxationFactor);

    std::string name;

private:
    SparseMatrix spMat_;
    mutable BiCGSTABIncompleteLUT solver_;
    SparseVector boundaries_, sources_;

    T& field_;
};

typedef Equation<ScalarFiniteVolumeField> ScalarEquation;
typedef Equation<VectorFiniteVolumeField> VectorEquation;

template<class T>
Equation<T> operator+(Equation<T> lhs, const Equation<T>& rhs);

template<class T>
Equation<T> operator-(Equation<T> lhs, const Equation<T>& rhs);

template<class T>
Equation<T> operator+(Equation<T> lhs, const T& rhs);

template<class T>
Equation<T> operator+(const T& lhs, Equation<T> rhs);

template<class T>
Equation<T> operator-(Equation<T> lhs, const T& rhs);

template<class T>
Equation<T> operator-(const T& lhs, Equation<T> rhs);

template<class T>
Equation<T> operator-(Equation<T> lhs, const T& rhs);

template<class T>
Equation<T> operator*(Equation<T> lhs, Scalar rhs);

template<class T>
Equation<T> operator*(Scalar lhs, Equation<T> rhs);

#include "Equation.tpp"

#endif
