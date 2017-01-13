#ifndef EQUATION_H
#define EQUATION_H

#include "FiniteVolumeGrid2D.h"
#include "ScalarFiniteVolumeField.h"
#include "VectorFiniteVolumeField.h"
#include "EigenSparseMatrixSolver.h"

template<class T>
class Equation
{
public:

    typedef SparseMatrixSolver::CoefficientList CoefficientList;

    //- Constructors
    Equation(T& field, const std::string& name = "N/A");
    Equation(const Input& input, T& field,
             const std::string& name,
             std::shared_ptr<SparseMatrixSolver>&& spSolver = std::make_shared<EigenSparseMatrixSolver>());
    Equation(const Equation<T>& other);

    //- Add/set/get coefficients
    void set(int i, int j, Scalar val);
    void add(int i, int j, Scalar val);

    Scalar& getRef(int i, int j);
    Scalar get(int i, int j) const;

    CoefficientList& coeffs() { return coeffs_; }
    const CoefficientList& coeffs() const { return coeffs_; }

    Vector& boundaries() { return boundaries_; }
    const Vector& boundaries() const { return boundaries_; }

    Vector& sources() { return sources_; }
    const Vector& sources() const { return sources_; }

    //- Clear equations
    void clear();

    //- Operators
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

    void setSparseSolver(std::shared_ptr<SparseMatrixSolver>& spSolver);
    std::shared_ptr<SparseMatrixSolver> & sparseSolver() { return spSolver_; }

    void configureSparseSolver(const Input& input);

    //- Solve the system
    Scalar solve();

    //- Relax the central coefficients (will fail if a central coefficient is not specified)
    void relax(Scalar relaxationFactor);

    std::string name;

private:

    CoefficientList coeffs_;
    Vector boundaries_, sources_;

    std::shared_ptr<SparseMatrixSolver> spSolver_;

    T& field_;
};

typedef Equation<ScalarFiniteVolumeField> ScalarEquation;
typedef Equation<VectorFiniteVolumeField> VectorEquation;

//- External functions
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
