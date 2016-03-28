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

    typedef Eigen::Triplet<Scalar> Triplet;

    Term() {}
    Term(const ScalarFiniteVolumeField& var);
    Term(const VectorFiniteVolumeField& var);

    Term& operator+=(const Term& rhs);
    Term& operator-=(const Term& rhs);
    Term& operator*=(Scalar rhs);
    Term& operator/=(Scalar rhs);
    Term& operator*=(const ScalarFiniteVolumeField& rhs);

    int nNonZeros() const { return coefficients_.size(); }
    const std::vector<Triplet>& coefficients() const { return coefficients_; }
    const std::vector<Scalar>& sources() const { return sources_; }

protected:

    std::vector<Triplet> coefficients_;
    std::vector<Scalar> sources_;

    friend Term operator==(Term term, Scalar rhs);
    friend Term operator==(Term term, const Term& rhs);
    friend Term operator==(Term term, const ScalarFiniteVolumeField& field);
};

Term operator==(Term term, Scalar rhs);
Term operator==(Term term, const Term& rhs);
Term operator==(Term term, const ScalarFiniteVolumeField& field);

Term operator+(Term lhs, const Term& rhs);
Term operator-(Term lhs, const Term& rhs);
Term operator*(Term lhs, Scalar rhs);
Term operator*(Scalar lhs, Term rhs);
Term operator*(const ScalarFiniteVolumeField& field, Term rhs);
Term operator*(Term lhs, const ScalarFiniteVolumeField& field);
Term operator/(Term lhs, Scalar rhs);

#endif
