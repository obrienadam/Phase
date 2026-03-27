#include <iostream>
#include <cmath>
#include <cassert>

#define TEST_ASSERT_EQ(a, b) \
    if ((a) != (b)) { \
        std::cerr << "FAIL: " << #a << " != " << #b << std::endl; \
        return 1; \
    }

#define TEST_ASSERT_NEAR(a, b, tol) \
    if (std::abs((a) - (b)) > (tol)) { \
        std::cerr << "FAIL: |" << #a << " - " << #b << "| > " << tol << std::endl; \
        return 1; \
    }

int main() {
    int tests_passed = 0;
    int tests_failed = 0;

    // Test 1: Vector length calculation
    {
        double x = 3.0, y = 4.0;
        double length = std::sqrt(x*x + y*y);
        if (length == 5.0) {
            std::cout << "[PASS] Vector length calculation" << std::endl;
            tests_passed++;
        } else {
            std::cout << "[FAIL] Vector length calculation" << std::endl;
            tests_failed++;
        }
    }

    // Test 2: Dot product
    {
        double ax = 1.0, ay = 2.0;
        double bx = 3.0, by = 4.0;
        double dotProduct = ax*bx + ay*by;
        if (dotProduct == 11.0) {
            std::cout << "[PASS] Dot product" << std::endl;
            tests_passed++;
        } else {
            std::cout << "[FAIL] Dot product" << std::endl;
            tests_failed++;
        }
    }

    // Test 3: Matrix determinant 2x2
    {
        double det = 1.0*4.0 - 2.0*3.0;
        if (det == -2.0) {
            std::cout << "[PASS] Matrix determinant 2x2" << std::endl;
            tests_passed++;
        } else {
            std::cout << "[FAIL] Matrix determinant 2x2" << std::endl;
            tests_failed++;
        }
    }

    // Test 4: Floating point tolerance
    {
        double a = 0.1 + 0.2;
        double b = 0.3;
        if (std::abs(a - b) < 1e-10) {
            std::cout << "[PASS] Floating point tolerance" << std::endl;
            tests_passed++;
        } else {
            std::cout << "[FAIL] Floating point tolerance" << std::endl;
            tests_failed++;
        }
    }

    // Test 5: Distance between points
    {
        double dist = std::sqrt(3.0*3.0 + 4.0*4.0);
        if (dist == 5.0) {
            std::cout << "[PASS] Distance between points" << std::endl;
            tests_passed++;
        } else {
            std::cout << "[FAIL] Distance between points" << std::endl;
            tests_failed++;
        }
    }

    std::cout << "\nMath Tests Summary: " << tests_passed << " passed, " << tests_failed << " failed" << std::endl;
    return tests_failed == 0 ? 0 : 1;
}

