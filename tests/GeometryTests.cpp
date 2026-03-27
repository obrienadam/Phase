#include <iostream>
#include <cmath>

int main() {
    int tests_passed = 0;
    int tests_failed = 0;

    // Test 1: 2D Vector magnitude
    {
        double x = 3.0, y = 4.0;
        double mag = std::sqrt(x*x + y*y);
        if (mag == 5.0) {
            std::cout << "[PASS] 2D Vector magnitude" << std::endl;
            tests_passed++;
        } else {
            std::cout << "[FAIL] 2D Vector magnitude" << std::endl;
            tests_failed++;
        }
    }

    // Test 2: Normalize 2D vector
    {
        double x = 3.0, y = 4.0;
        double mag = std::sqrt(x*x + y*y);
        double nx = x / mag;
        double ny = y / mag;
        double normalized_mag = std::sqrt(nx*nx + ny*ny);
        
        if (std::abs(normalized_mag - 1.0) < 1e-10) {
            std::cout << "[PASS] Normalize 2D vector" << std::endl;
            tests_passed++;
        } else {
            std::cout << "[FAIL] Normalize 2D vector" << std::endl;
            tests_failed++;
        }
    }

    // Test 3: Dot product 2D
    {
        double ax = 1.0, ay = 0.0;
        double bx = 0.0, by = 1.0;
        double dot = ax*bx + ay*by;
        
        if (dot == 0.0) {
            std::cout << "[PASS] Dot product 2D" << std::endl;
            tests_passed++;
        } else {
            std::cout << "[FAIL] Dot product 2D" << std::endl;
            tests_failed++;
        }
    }

    // Test 4: Cross product magnitude 2D
    {
        double ax = 1.0, ay = 0.0;
        double bx = 0.0, by = 1.0;
        double cross = ax*by - ay*bx;
        
        if (cross == 1.0) {
            std::cout << "[PASS] Cross product magnitude 2D" << std::endl;
            tests_passed++;
        } else {
            std::cout << "[FAIL] Cross product magnitude 2D" << std::endl;
            tests_failed++;
        }
    }

    // Test 5: Point in circle
    {
        double x = 3.0, y = 4.0;
        double dist = std::sqrt(x*x + y*y);
        double radius = 5.0;
        
        if (dist <= radius) {
            std::cout << "[PASS] Point in circle" << std::endl;
            tests_passed++;
        } else {
            std::cout << "[FAIL] Point in circle" << std::endl;
            tests_failed++;
        }
    }

    // Test 6: Bounding box area
    {
        double width = 10.0;
        double height = 5.0;
        double area = width * height;
        
        if (area == 50.0) {
            std::cout << "[PASS] Bounding box area" << std::endl;
            tests_passed++;
        } else {
            std::cout << "[FAIL] Bounding box area" << std::endl;
            tests_failed++;
        }
    }

    std::cout << "\nGeometry Tests Summary: " << tests_passed << " passed, " << tests_failed << " failed" << std::endl;
    return tests_failed == 0 ? 0 : 1;
}

