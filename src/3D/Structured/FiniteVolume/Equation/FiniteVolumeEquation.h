#ifndef PHASE_FINITE_VOLUME_EQUATION_H
#define PHASE_FINITE_VOLUME_EQUATION_H

#include "Math/Equation.h"
#include "Structured/FiniteVolume/Field/Field.h"
#include "Math/SparseMatrixSolver.h"

template<class T>
class FiniteVolumeEquation: public Equation
{
public:

    FiniteVolumeEquation(Field<T> &field);

    FiniteVolumeEquation(const std::string &name, Field<T> &field);

    FiniteVolumeEquation(const std::string &name, const Input &input, Field<T> &field);

    FiniteVolumeEquation& operator=(const Equation &rhs);

    void add(const Cell &cell, const Cell &nb, Scalar coeff);

    void addRhs(const Cell &cell, Scalar val);

    void configureSparseSolver(const Input &input);

    //- Solve

    Scalar solve() override;

protected:

        std::string _name;

        Field<T> &_field;
};

#include "FiniteVolumeEquation.tpp"

#endif
