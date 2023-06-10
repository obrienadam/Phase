#ifndef PHASE_FINITE_VOLUME_EQUATION_H
#define PHASE_FINITE_VOLUME_EQUATION_H

#include <valarray>

#include "System/Input.h"

#include "Math/CrsEquation.h"

#include "FiniteVolume/Field/ScalarFiniteVolumeField.h"
#include "FiniteVolume/Field/VectorFiniteVolumeField.h"
#include "FiniteVolumeGrid2D/FiniteVolumeGrid2D.h"

template <class T> class FiniteVolumeEquation : public CrsEquation {
public:
  //- Constructors

  FiniteVolumeEquation(FiniteVolumeField<T> &field,
                       const std::string &name = "", int nnz = 5);

  FiniteVolumeEquation(FiniteVolumeField<T> &field, int nnz)
      : FiniteVolumeEquation(field, "", nnz) {}

  FiniteVolumeEquation(const Input &input, FiniteVolumeField<T> &field,
                       const std::string &name, int nnz = 5);

  FiniteVolumeEquation(const FiniteVolumeEquation<T> &other) = default;

  FiniteVolumeEquation<T> &operator=(const FiniteVolumeEquation<T> &rhs);

  FiniteVolumeEquation<T> &operator=(FiniteVolumeEquation<T> &&rhs);

  FiniteVolumeEquation<T> &operator=(const CrsEquation &rhs);

  FiniteVolumeEquation<T> &operator=(CrsEquation &&rhs);

  //- Add/set/get coefficients
  void set(const Cell &cell, const Cell &nb, Scalar val);

  void add(const Cell &cell, const Cell &nb, Scalar val);

  void scale(const Cell &cell, Scalar val);

  template <class T2> void set(const Cell &cell, const Cell &nb, const T2 &val);

  template <class T2> void add(const Cell &cell, const Cell &nb, const T2 &val);

  template <typename cell_iterator, typename coeff_iterator>
  void add(const Cell &cell, cell_iterator begin, cell_iterator end,
           coeff_iterator coeffs) {
    for (cell_iterator itr = begin; itr != end; ++itr, ++coeffs)
      add(cell, *itr, *coeffs);
  }

  void add(const Cell &cell, const std::vector<Ref<const Cell>> &nbs,
           const std::vector<Scalar> &vals) {
    add(cell, nbs.begin(), nbs.end(), vals.begin());
  }

  void add(const Cell &cell, const std::vector<Ref<const Cell>> &nbs,
           const std::valarray<Scalar> &vals) {
    add(cell, nbs.begin(), nbs.end(), std::begin(vals));
  }

  T get(const Cell &cell, const Cell &nb);

  void remove(const Cell &cell);

  //- Set/get source vectors
  void addSource(const Cell &cell, Scalar val);

  void setSource(const Cell &cell, Scalar val);

  template <class T2> void addSource(const Cell &cell, const T2 &val);

  template <class T2> void setSource(const Cell &cell, const T2 &val);

  void configureSparseSolver(const Input &input, const Communicator &comm);

  //- Solve the system
  Scalar solve();

  //- Relax the central coefficients (will fail if a central coefficient is not
  //specified)
  void relax(Scalar relaxationFactor);

  const FiniteVolumeField<T> &field() const { return field_; }

  //- name
  std::string name;

protected:
  void mapFromSparseSolver();

  Size getRank() const;

  FiniteVolumeField<T> &field_;
};

template <class T>
FiniteVolumeEquation<T> operator*(const ScalarFiniteVolumeField &lhs,
                                  FiniteVolumeEquation<T> rhs);

#include "FiniteVolumeEquation.tpp"

#endif
