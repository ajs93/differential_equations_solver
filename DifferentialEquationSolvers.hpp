#ifndef __DIFFERENTIALECUATIONSOLVERS_HPP__
#define __DIFFERENTIALECUATIONSOLVERS_HPP__

#include <functional>
#include <utility>

#include <Matrix.hpp>

namespace DifferentialEquationSolvers {

template<typename T = double>
using InputGenerator = std::function<T(double t)>;

template<typename T = double>
using EcuationGenerator = std::function<T(double t, const T &y)>;

template<typename T = double>
using DiffPointSolver = std::function<T(double t_prev, const T &y_prev, double h, const EcuationGenerator<T> &diff_eq)>;

template<typename T = double>
struct DifferentialEcuationResult {
  std::vector<double> t;
  Matrix2D<T> y;

  DifferentialEcuationResult(const std::vector<double> &_t, const Matrix2D<T> &_y) : t(_t), y(_y) {}

  DifferentialEcuationResult(std::vector<double> &&_t, Matrix2D<T> &&_y) : t(std::move(_t)), y(std::move(_y)) {}

  DifferentialEcuationResult(const DifferentialEcuationResult &other) : t(other.t), y(other.y) {}

  DifferentialEcuationResult(DifferentialEcuationResult &&other) : t(std::move(other.t)), y(std::move(other.y)) {}
};

template<typename T = double>
T solveEulerPoint(double t_prev, const T &y_prev, double h, const EcuationGenerator<T> &diff_eq) {
  return y_prev + h * diff_eq(t_prev, y_prev);
}

template<typename T = double>
T solveHeunPoint(double t_prev, const T &y_prev, double h, const EcuationGenerator<T> &diff_eq) {
  T euler_term = diff_eq(t_prev, y_prev);
  T heun_term = diff_eq(t_prev + h, y_prev + h * diff_eq(t_prev, y_prev));
  return y_prev + ((h / 2.0) * (euler_term + heun_term));
}

template<typename T = double>
T solveRungeKuttaPoint(double t_prev, const T &y_prev, double h, const EcuationGenerator<T> &diff_eq) {
  auto f1 = diff_eq(t_prev, y_prev);
  auto f2 = diff_eq(t_prev + (h / 2.0), y_prev + ((h / 2.0) * f1));
  auto f3 = diff_eq(t_prev + (h / 2.0), y_prev + ((h / 2.0) * f2));
  auto f4 = diff_eq(t_prev + h, y_prev + (h * f3));

  return y_prev + h * ((f1 + (2 * f2) + (2 * f3) + f4) / 6);
}

template<typename T = double>
DifferentialEcuationResult<T> solve(const std::pair<double, double> &span,
                                    uint64_t sample_amount,
                                    const std::vector<T> &initial_conditions,
                                    const std::vector<EcuationGenerator<T>> &generators,
                                    const DiffPointSolver<T> &point_solver) {
  if (generators.size() != initial_conditions.size()) {
    throw std::runtime_error("Different generators and initial conditions size");
  }
  
  const double h = (span.second - span.first) / double(sample_amount - 1);
  std::vector<double> t(sample_amount);
  Matrix2D<T> y(sample_amount, generators.size());

  t.at(0) = span.first;

  for (uint64_t gen_idx = 0; gen_idx < generators.size(); gen_idx++) {
    y.at(0, gen_idx) = initial_conditions.at(gen_idx);
  }

  for (uint64_t t_idx = 1; t_idx < sample_amount; t_idx++) {
    t.at(t_idx) = span.first + (h * double(t_idx));

    for (uint64_t gen_idx = 0; gen_idx < generators.size(); gen_idx++) {
      y.at(t_idx, gen_idx) = point_solver(t.at(t_idx - 1), y.at(t_idx - 1, gen_idx), h, generators.at(gen_idx));
    }
  }

  return DifferentialEcuationResult(t, y);
}

template<typename T = double>
struct NthOrderDifferentialEcuationResult {
  std::vector<double> t;
  std::vector<Matrix2D<T>> y;

  NthOrderDifferentialEcuationResult(const std::vector<double> &_t, const std::vector<Matrix2D<T>> &_y) : t(_t), y(_y) {}

  NthOrderDifferentialEcuationResult(std::vector<double> &&_t, std::vector<Matrix2D<T>> &&_y) : t(std::move(_t)), y(std::move(_y)) {}

  NthOrderDifferentialEcuationResult(const NthOrderDifferentialEcuationResult &other) : t(other.t), y(other.y) {}

  NthOrderDifferentialEcuationResult(NthOrderDifferentialEcuationResult &&other) : t(std::move(other.t)), y(std::move(other.y)) {}
};


template<typename T = double>
NthOrderDifferentialEcuationResult<T> solveNthOrder(std::pair<double, double> span,
                                                    uint64_t sample_amount,
                                                    const std::vector<T> &initial_conditions,
                                                    const InputGenerator<T> &input_generator,
                                                    const DiffPointSolver<T> &point_solver,
                                                    const Matrix2D<T> &A,
                                                    const Matrix2D<T> &B) {
  if (A.getRowAmount() != A.getColumnAmount()) {
    throw std::runtime_error("A matrix must be a square matrix");
  }

  if (B.getColumnAmount() != 1) {
    throw std::runtime_error("B matrix must be a column matrix");
  }

  if (A.getRowAmount() != B.getRowAmount()) {
    throw std::runtime_error("A and B matrices must have the same row amount");
  }

  if (initial_conditions.size() != B.getRowAmount()) {
    throw std::runtime_error("Initial conditions and B row amount must be equal");
  }

  const uint64_t eq_amount = A.getRowAmount();
  const double h = (span.second - span.first) / double(sample_amount - 1);
  std::vector<double> t(sample_amount);
  std::vector<Matrix2D<T>> y;

  // Initial conditions for result matrix
  for (uint64_t eq_idx = 0; eq_idx < eq_amount; eq_idx++) {
    y.emplace_back(Matrix2D(eq_amount, 1));
    y.at(0).at(eq_idx, 0) = initial_conditions.at(eq_idx);
  }

  // Now we need to solve the N differential equations of first order
  for (uint64_t t_idx = 1; t_idx < sample_amount; t_idx++) {
    t.at(t_idx) = span.first + (h * double(t_idx));
    y.emplace_back(Matrix2D(eq_amount, 1));

    for (uint64_t eq_idx = 0; eq_idx < eq_amount; eq_idx++) {
      y.at(t_idx).at(eq_idx, 0) = point_solver(t.at(t_idx - 1), y.at(t_idx - 1).at(eq_idx, 0), h, [&y, &A, &B, &t_idx, &eq_idx, &input_generator] (double t, const T &_y) {
        auto a_x = A.getRow(eq_idx) * _y;
        auto b_u = B.getRow(eq_idx) * input_generator(t);
        return a_x.at(0, 0) + b_u.at(0, 0);
      });
    }
  }

  return NthOrderDifferentialEcuationResult(t, y);
}

};

#endif // __DIFFERENTIALECUATIONSOLVERS_HPP__