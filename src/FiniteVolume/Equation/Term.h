#ifndef TERM_H
#define TERM_H

#include <vector>
#include <eigen3/Eigen/Sparse>

#include "Types.h"
#include "SparseVector.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"
#include "FiniteVolumeGrid2D.h"

class Term
{
public:

    Term() {}
    Term(const FiniteVolumeGrid2D& grid);

    Term& operator+=(const Term& rhs);
    Term& operator-=(const Term& rhs);
    Term& operator*=(Scalar rhs);
    Term& operator/=(Scalar rhs);

protected:

    class Triplet : public Eigen::Triplet<Scalar>
    {
    public:

        using Eigen::Triplet<Scalar>::Triplet;

        Triplet& operator+=(const Triplet& other);
        Triplet& operator-=(const Triplet& other);
        Triplet& operator*=(Scalar rhs);
        Triplet& operator/=(Scalar rhs);
    };

    std::vector<Term::Triplet> coefficients_;
    std::vector<Scalar> sources_;
};

Term operator+(Term lhs, const Term& rhs);
Term operator-(Term lhs, const Term& rhs);
Term operator*(Term lhs, Scalar rhs);
Term operator*(Scalar lhs, Term rhs);
Term operator/(Term lhs, Scalar rhs);

Term div(const ScalarFiniteVolumeField& coeff, const ScalarFiniteVolumeField& var);
Term div(const ScalarFiniteVolumeField& coeff, const VectorFiniteVolumeField& var);
Term grad(const ScalarFiniteVolumeField& coeff, const ScalarFiniteVolumeField& var);
Term grad(const ScalarFiniteVolumeField& coeff, const VectorFiniteVolumeField& var);

#endif
