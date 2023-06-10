#ifndef PHASE_CRS_EQUATION_H
#define PHASE_CRS_EQUATION_H

#include "SparseMatrixSolver.h"

class CrsEquation {
public:
  CrsEquation() : rowPtr_(1, 0) {}

  CrsEquation(Size nRows, Size nnz);

  CrsEquation(const CrsEquation &eqn) = default;

  CrsEquation(CrsEquation &&eqn) = default;

  //- Assignment operators
  CrsEquation &operator=(const CrsEquation &eqn);

  CrsEquation &operator=(CrsEquation &&eqn);

  void setRank(Size rank);

  void setRank(Size nRows, Size nCols);

  void clear();

  Size rank() const { return rowPtr_.size() - 1; }

  Size capacity(Size row) const { return rowPtr_[row + 1] - rowPtr_[row]; }

  Size expand(Size row, Size nnz);

  Size expand(Size nnz);

  void addRow(Size nnz);

  void addRows(Size nRows, Size nnz);

  //- Add/set
  void addCoeff(Index localRow, Index globalCol, Scalar val);

  void setCoeff(Index localRow, Index globalCol, Scalar val);

  void scaleRow(Index localRow, Scalar val);

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
  const std::vector<Index> &rowPtr() const { return rowPtr_; }

  const std::vector<Index> &colInd() const { return colInd_; }

  const std::vector<Scalar> &vals() const { return vals_; }

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
  CrsEquation &operator+=(const CrsEquation &rhs);

  CrsEquation &operator-=(const CrsEquation &rhs);

  CrsEquation &operator+=(const Vector &rhs);

  CrsEquation &operator-=(const Vector &rhs);

  CrsEquation &operator*=(Scalar rhs);

  CrsEquation &operator/=(Scalar rhs);

  CrsEquation &operator==(Scalar rhs);

  CrsEquation &operator==(const CrsEquation &rhs);

  CrsEquation &operator==(const Vector &rhs);

protected:
  static std::vector<Index> tmpRowPtr_, tmpColInd_;

  static std::vector<Scalar> tmpVals_;

  std::vector<Index> rowPtr_, colInd_;

  std::vector<Scalar> vals_;

  Vector rhs_;

  std::shared_ptr<SparseMatrixSolver> solver_;
};

//- Operators

std::ostream &operator<<(std::ostream &os, const CrsEquation &eqn);

CrsEquation operator+(CrsEquation lhs, const CrsEquation &rhs);

CrsEquation operator-(CrsEquation lhs, const CrsEquation &rhs);

CrsEquation operator+(CrsEquation lhs, const Vector &rhs);

CrsEquation operator-(CrsEquation lhs, const Vector &rhs);

CrsEquation operator*(CrsEquation lhs, Scalar rhs);

CrsEquation operator/(CrsEquation lhs, Scalar rhs);

#endif
