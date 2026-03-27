#ifndef SIMPLE_TEST_FRAMEWORK_H
#define SIMPLE_TEST_FRAMEWORK_H

#include <iostream>
#include <vector>
#include <string>
#include <functional>
#include <cmath>

class SimpleTestFramework {
public:
    struct TestResult {
        std::string testName;
        bool passed;
        std::string message;
    };

    static SimpleTestFramework& instance() {
        static SimpleTestFramework framework;
        return framework;
    }

    void addTest(const std::string& name, std::function<bool()> testFunc) {
        tests.push_back({name, testFunc});
    }

    int runAll() {
        std::cout << "\n====== Running Unit Tests ======\n\n";
        int passed = 0;
        int failed = 0;

        for (const auto& test : tests) {
            try {
                bool result = test.second();
                if (result) {
                    std::cout << "[PASS] " << test.first << std::endl;
                    passed++;
                } else {
                    std::cout << "[FAIL] " << test.first << std::endl;
                    failed++;
                }
            } catch (const std::exception& e) {
                std::cout << "[ERROR] " << test.first << ": " << e.what() << std::endl;
                failed++;
            }
        }

        std::cout << "\n====== Test Summary ======\n";
        std::cout << "Passed: " << passed << "\n";
        std::cout << "Failed: " << failed << "\n";
        std::cout << "Total:  " << passed + failed << "\n\n";

        return failed == 0 ? 0 : 1;
    }

private:
    std::vector<std::pair<std::string, std::function<bool()>>> tests;
};

#define TEST_ASSERT_EQ(a, b) \
    if ((a) != (b)) { \
        std::cerr << "  Assertion failed: " << #a << " != " << #b << std::endl; \
        return false; \
    }

#define TEST_ASSERT_NEAR(a, b, tolerance) \
    if (std::abs((a) - (b)) > (tolerance)) { \
        std::cerr << "  Assertion failed: |" << #a << " - " << #b << "| > " << tolerance << std::endl; \
        return false; \
    }

#define TEST_ASSERT_TRUE(condition) \
    if (!(condition)) { \
        std::cerr << "  Assertion failed: " << #condition << " is false" << std::endl; \
        return false; \
    }

#define REGISTER_TEST(name, test_func) \
    namespace { \
        struct Registrar_##test_func { \
            Registrar_##test_func() { \
                SimpleTestFramework::instance().addTest(name, test_func); \
            } \
        } registrar_##test_func; \
    }

#endif // SIMPLE_TEST_FRAMEWORK_H
