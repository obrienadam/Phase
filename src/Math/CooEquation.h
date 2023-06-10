#ifndef PHASE_COO_EQUATION_H
#define PHASE_COO_EQUATION_H

#include "CrsEquation.h"
#include "SparseEntry.h"
#include "SparseMatrixSolver.h"

class CooEquation {
public:
  CooEquation() {}

  CooEquation(Size nRows, Size nnz);

  CooEquation(const CooEquation &eqn) = default;

  CooEquation(CooEquation &&eqn) = default;

  //- Assignment operators
  CooEquation &operator=(const CooEquation &eqn);

  CooEquation &operator=(CooEquation &&eqn);

  void setRank(Size rank);

  void setRank(Size nRows, Size nCols);

  Size rank() const { return rank_; }

  void clear();

  Size expand(Size row, Size nnz);

  Size expand(Size nnz);

  void addRow(Size nnz);

  void addRows(Size nRows, Size nnz);

  //- Add/set
  void addCoeff(Index localRow, Index globalCol, Scalar val);

  void setCoeff(Index localRow, Index globalCol, Scalar val);

  void addRhs(Index localRow, Scalar val) { rhs_(localRow) += val; }

  void setRhs(Index localRow, Scalar val) { rhs_(localRow) = val; }

  // - Multiple add/set

  template <class ColIt, class ValIt>
  void setCoeffs(Index row, ColIt colBegin, ColIt colEnd, ValIt valBegin) {
    while (colBegin != colEnd)
      setCoeff(row, *(colBegin++), *(valBegin++));
  }

  template <class ColIt, class ValIt>
  void addCoeffs(Index row, ColIt colBegin, ColIt colEnd, ValIt valBegin) {
    while (colBegin != colEnd)
      addCoeff(row, *(colBegin++), *(valBegin++));
  }

  void setCoeffs(Index row, const std::initializer_list<Index> &cols,
                 const std::initializer_list<Scalar> &vals) {
    setCoeffs(row, cols.begin(), cols.end(), vals.begin());
  }

  void addCoeffs(Index row, const std::initializer_list<Index> &cols,
                 const std::initializer_list<Scalar> &vals) {
    addCoeffs(row, cols.begin(), cols.end(), vals.begin());
  }

  //- Retrieve
  const std::vector<SparseEntry> &entries() const { return entries_; }

  Scalar coeff(Index localRow, Index globalCol) const;

  Scalar x(Index idx) const { return solver_->x(idx); }

  Scalar b(Index idx) const { return rhs_(idx); }

  //- Sparse solver

  void setSparseSolver(const std::shared_ptr<SparseMatrixSolver> &solver) {
    solver_ = solver;
  }

  const std::shared_ptr<SparseMatrixSolver> &sparseSolver() const {
    return solver_;
  }

  virtual Scalar solve();

  virtual Scalar solveLeastSquares();

  //- Operators

  CooEquation &operator+=(const std::vector<SparseEntry> &entries);

  CooEquation &operator+=(const CooEquation &rhs);

  CooEquation &operator-=(const CooEquation &rhs);

  CooEquation &operator+=(const Vector &rhs);

  CooEquation &operator-=(const Vector &rhs);

  CooEquation &operator*=(Scalar rhs);

  CooEquation &operator==(Scalar rhs);

  CooEquation &operator==(const CooEquation &rhs);

  CooEquation &operator==(const Vector &rhs);

protected:
  static std::vector<SparseEntry> tmp_;

  Size rank_;

  std::vector<SparseEntry> entries_;

  Vector rhs_;

  std::shared_ptr<SparseMatrixSolver> solver_;
};

//- Operators

std::ostream &operator<<(std::ostream &os, const CooEquation &eqn);

CooEquation operator+(CooEquation lhs, const CooEquation &rhs);

CooEquation operator-(CooEquation lhs, const CooEquation &rhs);

CooEquation operator+(CooEquation lhs, const Vector &rhs);

CooEquation operator-(CooEquation lhs, const Vector &rhs);

#endif
