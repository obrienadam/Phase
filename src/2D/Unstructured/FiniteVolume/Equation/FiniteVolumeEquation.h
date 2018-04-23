#ifndef PHASE_FINITE_VOLUME_EQUATION_H
#define PHASE_FINITE_VOLUME_EQUATION_H

#include <valarray>

#include "System/Input.h"

#include "Math/Equation.h"

#include "FiniteVolumeGrid2D/FiniteVolumeGrid2D.h"
#include "FiniteVolume/Field/ScalarFiniteVolumeField.h"
#include "FiniteVolume/Field/VectorFiniteVolumeField.h"

template<class T>
class FiniteVolumeEquation : public Equation
{
public:

    //- Constructors

    FiniteVolumeEquation(FiniteVolumeField<T> &field,
                         const std::string &name = "");

    FiniteVolumeEquation(const Input &input,
                         FiniteVolumeField<T> &field,
                         const std::string &name);

    Equation &operator =(const Equation &rhs) //- Necessary
    { return Equation::operator =(rhs); }

    Equation &operator =(Equation &&rhs)
    { return Equation::operator =(rhs); }

    //- Add/set/get coefficients
    void set(const Cell &cell, const Cell &nb, Scalar val);

    void add(const Cell &cell, const Cell &nb, Scalar val);

    template<class T2>
    void set(const Cell &cell, const Cell &nb, const T2 &val);

    template<class T2>
    void add(const Cell &cell, const Cell &nb, const T2 &val);

    template<typename cell_iterator, typename coeff_iterator>
    void add(const Cell &cell, cell_iterator begin, cell_iterator end, coeff_iterator coeffs)
    {
        for (cell_iterator itr = begin; itr != end; ++itr, ++coeffs)
            add(cell, *itr, *coeffs);
    }

    void add(const Cell &cell, const std::vector<Ref<const Cell>> &nbs, const std::vector<Scalar> &vals)
    {add(cell, nbs.begin(), nbs.end(), vals.begin());}

    void add(const Cell &cell, const std::vector<Ref<const Cell>> &nbs, const std::valarray<Scalar> &vals)
    {add(cell, nbs.begin(), nbs.end(), std::begin(vals));}

    void addCoupling(const Cell &cell, const Cell &nb, const T &val);

    T get(const Cell &cell, const Cell &nb);

    void remove(const Cell &cell);

    //- Set/get source vectors
    void addSource(const Cell &cell, Scalar val);

    void setSource(const Cell &cell, Scalar val);

    template<class T2>
    void addSource(const Cell &cell, const T2 &val);

    template<class T2>
    void setSource(const Cell &cell, const T2 &val);

    void setSparseSolver(const std::shared_ptr<SparseMatrixSolver> &spSolver);

    void configureSparseSolver(const Input &input, const Communicator &comm);

    //- Solve the system
    Scalar solve();

    //- Relax the central coefficients (will fail if a central coefficient is not specified)
    void relax(Scalar relaxationFactor);

    //- name
    std::string name;

protected:

    void mapFromSparseSolver();

    Size getRank() const;

    FiniteVolumeField<T> &field_;
};

#include "FiniteVolumeEquation.tpp"

#endif
