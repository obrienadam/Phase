#define BOOST_TEST_MODULE MatrixTest
#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include "Matrix.h"
#include "ImmersedBoundary/BilinearInterpolation.h"

BOOST_AUTO_TEST_SUITE (MatrixTest)

BOOST_AUTO_TEST_CASE (AllocationTest)
{
    Matrix A = random(10, 8, 0, 15);
    Matrix B = A;

    BOOST_REQUIRE_EQUAL(A.nRows(), 10);
    BOOST_REQUIRE_EQUAL(A.nCols(), 8);
    BOOST_REQUIRE_EQUAL(A.nRows(), B.nRows());
    BOOST_REQUIRE_EQUAL(A.nCols(), B.nCols());
    BOOST_REQUIRE(A.data() != B.data());

    for(int i = 0; i < A.nRows(); ++i)
        for(int j = 0; j < A.nCols(); ++j)
            BOOST_REQUIRE_EQUAL(A(i, j), B(i, j));
}

BOOST_AUTO_TEST_CASE (TransposeTest)
{
    // Square matrix
    Matrix A = random(15, 15, 0, 15);
    Matrix B = transpose(A);

    for(int i = 0; i < A.nRows(); ++i)
        for(int j = 0; j < A.nCols(); ++j)
            BOOST_REQUIRE_EQUAL(A(i, j), B(j, i));

    // Non-square matrix
    A = random(15, 10, 0, 15);
    B = transpose(A);

    for(int i = 0; i < A.nRows(); ++i)
        for(int j = 0; j < A.nCols(); ++j)
            BOOST_REQUIRE_EQUAL(A(i, j), B(j, i));
}

BOOST_AUTO_TEST_CASE (SolverTest)
{
    Matrix A = random(15, 15, 0, 10);
    Matrix B = random(15, 2, 0, 10);
    Matrix X1 = solve(A, B);
    Matrix X2 = inverse(A)*B;
    Matrix I = inverse(A)*A;

    for(int i = 0; i < X1.nRows(); ++i)
        for(int j = 0; j < X1.nCols(); ++j)
            BOOST_REQUIRE_CLOSE(X1(i, j), X2(i, j), 1e-8);

    for(int i = 0; i < I.nRows(); ++i)
        for(int j = 0; j < I.nCols(); ++j)
        {
            if(i == j)
                BOOST_REQUIRE_CLOSE(I(i, j), 1., 1e-9);
            else
                BOOST_REQUIRE_CLOSE(I(i, j) + 1., 1., 1e-9);
        }
}

BOOST_AUTO_TEST_CASE(LeastSquaresTest)
{
    Matrix A = random(18, 10, 0, 10);
    Matrix B = random(18, 2, 0, 10);
    Matrix X1 = solve(A, B);
    Matrix X2 = inverse(transpose(A)*A)*transpose(A)*B;

    for(int i = 0; i < X1.nRows(); ++i)
        for(int j = 0; j < X1.nCols(); ++j)
            BOOST_REQUIRE_CLOSE(X1(i, j), X2(i, j), 1e-9);
}

BOOST_AUTO_TEST_CASE(BilinearInterpolationTest)
{
    Point2D pts[] = {
        Point2D(0, 0),
        Point2D(1.2, 0),
        Point2D(0.9, 1.3),
        Point2D(0, 0.86),
    };

    Point2D ip(0.5, 1.1);

    Scalar vals[] = {0.1, 1, 1, 0.5};

    BilinearInterpolation bl(pts);

    auto coeffs = bl(ip);

    Scalar sumCoeff = 0.;

    for(auto coeff: coeffs)
        sumCoeff += coeff;

    BOOST_REQUIRE_CLOSE(1., sumCoeff, 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()
