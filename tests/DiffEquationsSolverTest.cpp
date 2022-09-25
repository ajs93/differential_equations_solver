#include <catch2/catch_test_macros.hpp>

#include <catch2/catch_approx.hpp>

#include <DifferentialEquationSolvers.hpp>
#include <math.h>

SCENARIO("Simple test example", "[auto]") {
  // Example taken from: https://en.wikipedia.org/wiki/Euler_method
  GIVEN("A differential ecuation with known solution, an initial condition and a time span") {
    DifferentialEquationSolvers::EcuationGenerator<double> generator = [](double t, const double &y) {
      return y;
    };

    double initial_condition = 1.0;
    std::pair<double, double> span = {0, 4};
    uint64_t sample_amount = 5;

    WHEN("The differential ecuation solver is executed with known values") {
      auto res = DifferentialEquationSolvers::solve<double>(span, sample_amount, {initial_condition}, {generator}, DifferentialEquationSolvers::solveEulerPoint<double>);

      THEN("The results are well known") {
        REQUIRE(res.t.at(0) == Catch::Approx(0).epsilon(0.01));
        REQUIRE(res.t.at(1) == Catch::Approx(1).epsilon(0.01));
        REQUIRE(res.t.at(2) == Catch::Approx(2).epsilon(0.01));
        REQUIRE(res.t.at(3) == Catch::Approx(3).epsilon(0.01));
        REQUIRE(res.t.at(4) == Catch::Approx(4).epsilon(0.01));

        REQUIRE(res.y.at(0).at(0) == initial_condition);
        REQUIRE(res.y.at(0).at(1) == Catch::Approx(2).epsilon(0.01));
        REQUIRE(res.y.at(0).at(2) == Catch::Approx(4).epsilon(0.01));
        REQUIRE(res.y.at(0).at(3) == Catch::Approx(8).epsilon(0.01));
        REQUIRE(res.y.at(0).at(4) == Catch::Approx(16).epsilon(0.01));
      }
    }
  }

  // Example taken from: https://tutorial.math.lamar.edu/classes/de/eulersmethod.aspx
  GIVEN("A differential ecuation (IVP) with known solution, an initial condition and a time span") {
    DifferentialEquationSolvers::EcuationGenerator<double> ivp_gen = [](double t, const double &y) {
      return 2 - exp((-4) * t) - (2 * y);
    };

    double initial_condition = 1.0;
    std::pair<double, double> span = {0, 0.5};
    uint64_t sample_amount = 6;

    WHEN("The differential ecuation solver is executed with known values") {
      auto res = DifferentialEquationSolvers::solve<double>(span, sample_amount, {initial_condition}, {ivp_gen}, DifferentialEquationSolvers::solveEulerPoint<double>);

      THEN("The results are well known") {
        REQUIRE(res.t.at(0) == Catch::Approx(0).epsilon(0.01));
        REQUIRE(res.t.at(1) == Catch::Approx(0.1).epsilon(0.01));
        REQUIRE(res.t.at(2) == Catch::Approx(0.2).epsilon(0.01));
        REQUIRE(res.t.at(3) == Catch::Approx(0.3).epsilon(0.01));
        REQUIRE(res.t.at(4) == Catch::Approx(0.4).epsilon(0.01));
        REQUIRE(res.t.at(5) == Catch::Approx(0.5).epsilon(0.01));

        REQUIRE(res.y.at(0).at(0) == Catch::Approx(1).epsilon(0.01));
        REQUIRE(res.y.at(0).at(1) == Catch::Approx(0.9).epsilon(0.01));
        REQUIRE(res.y.at(0).at(2) == Catch::Approx(0.852967995).epsilon(0.01));
        REQUIRE(res.y.at(0).at(3) == Catch::Approx(0.837441500).epsilon(0.01));
        REQUIRE(res.y.at(0).at(4) == Catch::Approx(0.839833779).epsilon(0.01));
        REQUIRE(res.y.at(0).at(5) == Catch::Approx(0.851677371).epsilon(0.01));
      }
    }
  }
}

SCENARIO("N-th order differential equation with 1x1 matrix", "[auto]") {
  // Reutilizing IVP example that we solved using euler equation
  GIVEN("A differential equation system of order 1 with well known values") {
    DifferentialEquationSolvers::InputGenerator<double> input_generator = [](double t) {
      return 2 - exp((-4) * t);
    };

    std::vector<double> initial_conditions = {1.0};
    std::pair<double, double> span = {0, 0.5};
    uint64_t sample_amount = 6;

    Matrix2D A(std::vector<std::vector<double>>{{-2}});
    Matrix2D B(std::vector<std::vector<double>>{{1}}); 

    WHEN("The differential equations solver is executed with euler method") {
      auto res = DifferentialEquationSolvers::solveNthOrder<double>(span, sample_amount, initial_conditions, input_generator, DifferentialEquationSolvers::solveEulerPoint<double>, A, B);

      THEN("Values check") {
        REQUIRE(res.t.at(0) == Catch::Approx(0).epsilon(0.01));
        REQUIRE(res.t.at(1) == Catch::Approx(0.1).epsilon(0.01));
        REQUIRE(res.t.at(2) == Catch::Approx(0.2).epsilon(0.01));
        REQUIRE(res.t.at(3) == Catch::Approx(0.3).epsilon(0.01));
        REQUIRE(res.t.at(4) == Catch::Approx(0.4).epsilon(0.01));
        REQUIRE(res.t.at(5) == Catch::Approx(0.5).epsilon(0.01));

        REQUIRE(res.y.at(0).at(0) == Catch::Approx(1).epsilon(0.01));
        REQUIRE(res.y.at(0).at(1) == Catch::Approx(0.9).epsilon(0.01));
        REQUIRE(res.y.at(0).at(2) == Catch::Approx(0.852967995).epsilon(0.01));
        REQUIRE(res.y.at(0).at(3) == Catch::Approx(0.837441500).epsilon(0.01));
        REQUIRE(res.y.at(0).at(4) == Catch::Approx(0.839833779).epsilon(0.01));
        REQUIRE(res.y.at(0).at(5) == Catch::Approx(0.851677371).epsilon(0.01));
      }
    }
  }
}