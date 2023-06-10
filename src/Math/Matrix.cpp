#include <algorithm>
#include <ostream>

#ifdef __INTEL_COMPILER
#include <mkl.h>
#else
extern "C" {
#define lapack_complex_float                                                   \
  float _Complex // Fixes the world's most ridiculous bug
#define lapack_complex_double double _Complex
#include <cblas.h>
#include <lapacke.h>
}
#endif

#include "System/Exception.h"

#include "Matrix.h"

Matrix Matrix::_tmp;

Matrix::Matrix(Size m, Size n, const std::initializer_list<Scalar> &coeffs) {
  vals_.assign(coeffs);
  resize(m, n);
}

void Matrix::resize(Size m, Size n) {
  m_ = m;
  n_ = n;
  vals_.resize(m_ * n_, 0.);
  ipiv_.resize(m_);
}

void Matrix::zero() { std::fill(vals_.begin(), vals_.end(), 0.); }

void Matrix::setIdentity(Size m) {
  resize(m, m);
  for (auto i = 0; i < m_; ++i)
    for (auto j = 0; j < n_; ++j)
      vals_[i * n_ + j] = i == j ? 1. : 0.;
}

void Matrix::init(const Scalar *begin, const Scalar *end) {
  std::copy(begin, end, vals_.begin());
}

void Matrix::setRow(int i, const std::initializer_list<Scalar> &coeffs) {
  std::copy(coeffs.begin(), coeffs.end(), vals_.begin() + n_ * i);
}

void Matrix::addToRow(int i, const std::initializer_list<Scalar> &coeffs) {
  std::transform(coeffs.begin(), coeffs.end(), vals_.begin() + n_ * i,
                 vals_.begin() + n_ * i, std::plus<Scalar>());
}

Scalar &Matrix::operator()(Size i, Size j) { return vals_[i * n_ + j]; }

Scalar Matrix::operator()(Size i, Size j) const { return vals_[i * n_ + j]; }

Matrix &Matrix::operator=(const std::initializer_list<Scalar> &coeffs) {
  if (m_ * n_ != coeffs.size())
    throw Exception("Matrix", "operator=", "dimension mismatch.");

  vals_.assign(coeffs);

  return *this;
}

//- Operators
Matrix &Matrix::operator+=(const Matrix &rhs) {
  if (m_ != rhs.m_ || n_ != rhs.n_)
    throw Exception("Matrix", "operator+=", "dimension mismatch.");

  std::transform(vals_.begin(), vals_.end(), rhs.vals_.begin(), vals_.begin(),
                 std::plus<Scalar>());

  return *this;
}

Matrix &Matrix::operator-=(const Matrix &rhs) {
  if (m_ != rhs.m_ || n_ != rhs.n_)
    throw Exception("Matrix", "operator-=", "dimension mismatch.");

  std::transform(vals_.begin(), vals_.end(), rhs.vals_.begin(), vals_.begin(),
                 std::minus<Scalar>());

  return *this;
}

Matrix &Matrix::operator*=(Scalar rhs) {
  std::for_each(vals_.begin(), vals_.end(), [rhs](Scalar &val) { val *= rhs; });
  return *this;
}

Matrix &Matrix::operator/=(Scalar rhs) {
  std::for_each(vals_.begin(), vals_.end(), [rhs](Scalar &val) { val /= rhs; });
  return *this;
}

Matrix &Matrix::solve(Matrix &b) {
  if (isSquare())
    LAPACKE_dgesv(LAPACK_ROW_MAJOR, m_, b.n_, data(), n_, ipiv_.data(),
                  b.data(), b.n_);
  else {
    LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', m_, n_, b.n_, data(), n_, b.data(),
                  b.n_);
    b.m_ = n_; // A bit hackish, but resizes b apropriately. Only works so long
               // as the data format is row major
  }

  return b;
}

Matrix Matrix::solve(const Matrix &b) {
  Matrix x(b.m(), b.n());
  solve(x);

  return x;
}

Matrix &Matrix::transpose() {
  if (isSquare()) // Square matrices
  {
    auto &self = *this;

    for (size_t m = 0; m < m_; ++m)
      for (size_t n = m + 1; n < n_; ++n)
        std::swap(self(m, n), self(n, m));
  } else // General non-square matrices
  {
    auto last = vals_.end();
    auto first = vals_.begin();
    int m = n_;

    const int mn1 = (last - first - 1);
    const int n = (last - first) / m;
    std::vector<bool> visited(last - first);
    auto cycle = first;
    while (++cycle != last) {
      if (visited[cycle - first])
        continue;
      int a = cycle - first;
      do {
        a = a == mn1 ? mn1 : (n * a) % mn1;
        std::swap(*(first + a), *cycle);
        visited[a] = true;
      } while ((first + a) != cycle);
    }

    std::swap(m_, n_);
  }

  return *this;
}

Matrix &Matrix::invert() {
  lapack_int info1 =
      LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m_, n_, data(), n_, ipiv_.data());
  lapack_int info2 =
      LAPACKE_dgetri(LAPACK_ROW_MAJOR, m_, data(), n_, ipiv_.data());

  if (info1 != 0 || info2 != 0)
    throw Exception(
        "Matrix", "invert",
        "inversion failed, matrix is singular to working precision.");

  return *this;
}

Matrix &Matrix::pinvert() {
  _tmp.setIdentity(std::max(m_, n_));
  LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', m_, n_, _tmp.n(), data(), n_,
                _tmp.data(), _tmp.n());
  _tmp.resize(n_, m_);
  return (*this = _tmp);
}

Scalar Matrix::norm(char type) const {
  return LAPACKE_dlange(LAPACK_ROW_MAJOR, type, m_, n_, data(), n_);
}

Scalar Matrix::cond(char type) const {
  Scalar rcond;
  LAPACKE_dgecon(LAPACK_ROW_MAJOR, type, m_, data(), n_, norm(type), &rcond);
  return rcond;
}

Matrix Matrix::subMatrix(size_t startRow, size_t startCol, size_t endRow,
                         size_t endCol) const {
  Matrix subMat(endRow - startRow, endCol - startCol);
  auto &self = *this;

  for (size_t m = startRow; m < endRow; ++m)
    for (size_t n = startCol; n < endCol; ++n)
      subMat(m - startRow, n - startCol) = self(m, n);

  return subMat;
}

//- External functions
Matrix eye(int m) {
  Matrix mat(m, m);
  for (int i = 0; i < m; ++i)
    mat(i, i) = 1.;

  return mat;
}

Matrix random(size_t nRows, size_t nCols, Scalar min, Scalar max) {
  Matrix mat(nRows, nCols);

  for (size_t m = 0; m < nRows; ++m)
    for (size_t n = 0; n < nCols; ++n)
      mat(m, n) = (double)rand() / RAND_MAX * (max - min) + min;

  return mat;
}

Matrix transpose(const Matrix &mat) {
  Matrix transMat(mat.n(), mat.m());

  for (int i = 0; i < mat.m(); ++i)
    for (int j = 0; j < mat.n(); ++j)
      transMat(j, i) = mat(i, j);

  return transMat;
}

Matrix inverse(Matrix mat) { return mat.invert(); }

Matrix pseudoInverse(Matrix mat) {
  Matrix I = eye(std::max(mat.m(), mat.n()));
  LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', mat.m(), mat.n(), I.n(), mat.data(),
                mat.n(), I.data(), I.n());
  I.resize(mat.n(), mat.m());
  return I;
}

Matrix solve(Matrix A, Matrix b) { return A.solve(b); }

Matrix operator+(Matrix lhs, const Matrix &rhs) { return lhs += rhs; }

Matrix operator-(Matrix lhs, const Matrix &rhs) { return lhs -= rhs; }

Matrix operator*(Matrix lhs, Scalar rhs) { return lhs *= rhs; }

Matrix operator*(Scalar lhs, Matrix rhs) { return rhs *= lhs; }

Matrix operator*(const Matrix &A, const Matrix &B) {
  Matrix C(A.m(), B.n());

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A.m(), B.n(), A.n(),
              1., A.data(), A.n(), B.data(), B.n(), 1., C.data(), C.n());

  // Works for sure but slower
  //    const int nI = A.nRows(), nJ = B.nCols(), nK = A.nCols();

  //    for(int i = 0; i < nI; ++i)
  //        for(int j = 0; j < nJ; ++j)
  //            for(int k = 0; k < nK; ++k)
  //                C(i, j) += A(i, k)*B(k, j);

  return C;
}

Matrix operator/(Matrix lhs, Scalar rhs) { return lhs /= rhs; }

Matrix multiply(const Matrix &A, const Matrix &B, bool transA, bool transB) {
  Matrix C(A.m(), B.n());

  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, A.m(), B.n(), A.n(), 1.,
              A.data(), A.n(), B.data(), B.n(), 1., C.data(), C.n());
  //
  //    cblas_dgemm(CblasRowMajor, transA ? CblasTrans : CblasNoTrans, transB ?
  //    CblasTrans : CblasNoTrans,
  //                A.m(), B.n(), A.n(), 1., A.data(), A.n(),
  //                B.data(), B.n(), 1., C.data(), C.n());

  return C;
}

std::ostream &operator<<(std::ostream &os, const Matrix &mat) {
  for (Size i = 0, m = mat.m(); i < m; ++i) {
    for (Size j = 0, n = mat.n(); j < n; ++j) {
      os << mat(i, j) << ' ';
    }
    os << '\n';
  }

  return os;
}
