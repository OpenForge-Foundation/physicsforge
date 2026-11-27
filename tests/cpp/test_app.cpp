// tests/cpp/test_app.cpp - entry point for the physicsforge C++ unit tests.
#include <iostream>

#include "test_framework.h"

void test_constants();
void test_vector();
void test_mechanics();
void test_electromagnetism();
void test_thermodynamics();
void test_quantum();
void test_relativity();
void test_numerics();
void test_ode();
void test_api();

int main() {
    std::cout << "Running physicsforge C++ tests...\n";

    struct Suite {
        const char *name;
        void (*run)();
    };
    constexpr Suite suites[] = {
            {"constants", test_constants},
            {"vector", test_vector},
            {"mechanics", test_mechanics},
            {"electromagnetism", test_electromagnetism},
            {"thermodynamics", test_thermodynamics},
            {"quantum", test_quantum},
            {"relativity", test_relativity},
            {"numerics", test_numerics},
            {"ode", test_ode},
            {"c api", test_api},
    };

    for (const auto &suite: suites) {
        const int before = pftest::failure_count();
        suite.run();
        const int failed = pftest::failure_count() - before;
        std::cout << "  [" << (failed == 0 ? "PASS" : "FAIL") << "] " << suite.name << "\n";
    }

    std::cout << pftest::check_count() << " checks, " << pftest::failure_count() << " failures\n";
    return pftest::failure_count() == 0 ? 0 : 1;
}
