#include <catch2/catch_test_macros.hpp>

#include <catch2/catch_approx.hpp>

#include <DifferentialEquationSolvers.hpp>
#include <math.h>

SCENARIO("Ejercicio 2.a", "[auto]") {
  GIVEN("A differential equation, initial condition, a span and sample amount with the real solution") {
    DifferentialEquationSolvers::EcuationGenerator<double> gen = [](double t, const double &y) {
      return pow(t, 2) - y;
    };

    DifferentialEquationSolvers::InputGenerator<double> real_solution = [](double t) {
      return -1 * exp(-1 * t) + pow(t, 2) - 2 * t + 2;
    };

    double initial_condition = 1.0;
    std::pair<double, double> span = {0, 2};
    uint64_t sample_amount = 20;

    WHEN("Euler method is used to solve the equation") {
      auto ret = DifferentialEquationSolvers::solve<double>(span, sample_amount, {initial_condition}, {gen}, DifferentialEquationSolvers::solveEulerPoint<double>);

      THEN("The error shall be inside 10%% of the real solution") {
        for (uint64_t t_idx = 0; t_idx < ret.t.size(); t_idx++) {
          CHECK(ret.y.at(t_idx, 0) == Catch::Approx(real_solution(ret.t.at(t_idx))).epsilon(0.1));
        }
      }
    }

    WHEN("Heun method is used to solve the equation") {
      auto ret = DifferentialEquationSolvers::solve<double>(span, sample_amount, {initial_condition}, {gen}, DifferentialEquationSolvers::solveHeunPoint<double>);

      THEN("The error shall be inside 1%% of the real solution") {
        for (uint64_t t_idx = 0; t_idx < ret.t.size(); t_idx++) {
          CHECK(ret.y.at(t_idx, 0) == Catch::Approx(real_solution(ret.t.at(t_idx))).epsilon(0.01));
        }
      }
    }

    WHEN("Runge kutta is used to solve the equation") {
      auto ret = DifferentialEquationSolvers::solve<double>(span, sample_amount, {initial_condition}, {gen}, DifferentialEquationSolvers::solveRungeKuttaPoint<double>);

      THEN("The error shall be inside 0.1%% of the real solution") {
        for (uint64_t t_idx = 0; t_idx < ret.t.size(); t_idx++) {
          CHECK(ret.y.at(t_idx, 0) == Catch::Approx(real_solution(ret.t.at(t_idx))).epsilon(0.001));
        }
      }
    }
  }
}