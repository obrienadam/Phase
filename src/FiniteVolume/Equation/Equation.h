#ifndef EQUATION_H
#define EQUATION_H

#include "FiniteVolumeGrid2D.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"
#include "SparseMatrixSolver.h"
#include "Communicator.h"

template<class T>
class Equation
{
public:

    typedef SparseMatrixSolver::CoefficientList CoefficientList;

    //- Constructors
    Equation(FiniteVolumeField<T> &field,
             const std::string &name = "N/A");

    Equation(const Input &input,
             FiniteVolumeField<T> &field,
             const std::string &name);

    Equation(const Equation<T> &rhs) = default;

    Equation(Equation<T> &&rhs) = default;

    //- Add/set/get coefficients
    template<typename T2>
    void set(const Cell &cell, const Cell &nb, T2 val);

    template<typename T2>
    void add(const Cell &cell, const Cell &nb, T2 val);

    T get(const Cell &cell, const Cell &nb);

    void remove(const Cell &cell);

    //- Get a coefficient adjacency list, used only for assembling systems
    const CoefficientList &coeffs() const
    { return coeffs_; }

    //- Set/get source vectors
    void addSource(const Cell &cell, T val);

    void setSource(const Cell &cell, T val);

    const Vector &sources() const
    { return sources_; }

    //- Clear equations
    void clear();

    //- Operators
    Equation<T> &operator=(const Equation<T> &rhs);

    Equation<T> &operator=(Equation<T> &&rhs);

    Equation<T> &operator+=(const Equation<T> &rhs);

    Equation<T> &operator-=(const Equation<T> &rhs);

    Equation<T> &operator+=(const FiniteVolumeField<T> &rhs);

    Equation<T> &operator-=(const FiniteVolumeField<T> &rhs);

    Equation<T> &operator*=(Scalar rhs);

    Equation<T> &operator*=(const ScalarFiniteVolumeField &rhs);

    Equation<T> &operator/=(const ScalarFiniteVolumeField& rhs);

    Equation<T> &operator==(Scalar rhs);

    Equation<T> &operator==(const Equation<T> &rhs);

    Equation<T> &operator==(const FiniteVolumeField<T> &rhs);

    void setSparseSolver(std::shared_ptr<SparseMatrixSolver> &spSolver);

    std::shared_ptr<SparseMatrixSolver> &sparseSolver()
    { return spSolver_; }

    void configureSparseSolver(const Input &input, const Communicator &comm);

    //- Solve the system
    Scalar solve();

    Scalar solveWithGuess();

    //- Relax the central coefficients (will fail if a central coefficient is not specified)
    void relax(Scalar relaxationFactor);

    //- name
    std::string name;

protected:

    Size getRank() const;

    void setValue(Index i, Index j, Scalar val);

    void addValue(Index i, Index j, Scalar val);

    Scalar &coeffRef(Index i, Index j);

    Size nActiveCells_; // Cached for efficiency

    CoefficientList coeffs_;
    Vector sources_;

    std::shared_ptr<SparseMatrixSolver> spSolver_;

    FiniteVolumeField<T> &field_;
};

//- External functions
template<class T>
Equation<T> operator+(Equation<T> lhs, const Equation<T> &rhs);

template<class T>
Equation<T> operator-(Equation<T> lhs, const Equation<T> &rhs);

template<class T>
Equation<T> operator+(Equation<T> lhs, const FiniteVolumeField<T> &rhs);

template<class T>
Equation<T> operator+(const FiniteVolumeField<T> &lhs, Equation<T> rhs);

template<class T>
Equation<T> operator-(Equation<T> lhs, const FiniteVolumeField<T> &rhs);

template<class T>
Equation<T> operator-(const FiniteVolumeField<T> &lhs, Equation<T> rhs);

template<class T>
Equation<T> operator-(Equation<T> lhs, const FiniteVolumeField<T> &rhs);

template<class T>
Equation<T> operator*(Equation<T> lhs, Scalar rhs);

template<class T>
Equation<T> operator*(Scalar lhs, Equation<T> rhs);

template<class T>
Equation<T> operator/(Equation<T> lhs, const ScalarFiniteVolumeField& rhs)
{
    return lhs /= rhs;
}

#include "Equation.tpp"

#endif
