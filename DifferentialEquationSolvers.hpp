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

/**
 * @brief Result from resolving a single order differential equation
 * 
 * @tparam T Template parameter for the data type of the result
 */
template<typename T = double>
struct DifferentialEquationResult {
  /** Time vector for the result */
  std::vector<double> t;

  /**
   * @brief Matrix of results where each column represents
   * a result vector
   */
  Matrix2D<T> y;

  DifferentialEquationResult(const std::vector<double> &_t, const Matrix2D<T> &_y) : t(_t), y(_y) {}

  DifferentialEquationResult(std::vector<double> &&_t, Matrix2D<T> &&_y) : t(std::move(_t)), y(std::move(_y)) {}

  DifferentialEquationResult(const DifferentialEquationResult &other) : t(other.t), y(other.y) {}

  DifferentialEquationResult(DifferentialEquationResult &&other) : t(std::move(other.t)), y(std::move(other.y)) {}
};

/**
 * @brief Calculate a differential point with the Euler method
 * 
 * @tparam T Template data type for this solver
 * @param t_prev Previous time for the calculation
 * @param y_prev Previous result for the calculation
 * @param h Sample time delta
 * @param diff_eq Differential equation generator
 * @return T Calculated point
 */
template<typename T = double>
T solveEulerPoint(double t_prev, const T &y_prev, double h, const EcuationGenerator<T> &diff_eq) {
  return y_prev + h * diff_eq(t_prev, y_prev);
}

/**
 * @brief Calculate a differential point with the Heun/Newton method
 * 
 * @tparam T Template data type for this solver
 * @param t_prev Previous time for the calculation
 * @param y_prev Previous result for the calculation
 * @param h Sample time delta
 * @param diff_eq Differential equation generator
 * @return T Calculated point
 */
template<typename T = double>
T solveHeunPoint(double t_prev, const T &y_prev, double h, const EcuationGenerator<T> &diff_eq) {
  T euler_term = diff_eq(t_prev, y_prev);
  T heun_term = diff_eq(t_prev + h, y_prev + h * diff_eq(t_prev, y_prev));
  return y_prev + ((h / 2.0) * (euler_term + heun_term));
}

/**
 * @brief Calculate a differential point with the Runge-Kutta 4 method
 * 
 * @tparam T Template data type for this solver
 * @param t_prev Previous time for the calculation
 * @param y_prev Previous result for the calculation
 * @param h Sample time delta
 * @param diff_eq Differential equation generator
 * @return T Calculated point
 */
template<typename T = double>
T solveRungeKuttaPoint(double t_prev, const T &y_prev, double h, const EcuationGenerator<T> &diff_eq) {
  auto f1 = diff_eq(t_prev, y_prev);
  auto f2 = diff_eq(t_prev + (h / 2.0), y_prev + ((h / 2.0) * f1));
  auto f3 = diff_eq(t_prev + (h / 2.0), y_prev + ((h / 2.0) * f2));
  auto f4 = diff_eq(t_prev + h, y_prev + (h * f3));

  return y_prev + h * ((f1 + (2 * f2) + (2 * f3) + f4) / 6);
}

/**
 * @brief Solve multiple independent single order differential equations
 * @exception std::runtime_error If there is a mismatch between generators and initial
 * conditions vector sizes
 * 
 * @tparam T Template data type for this solver
 * @param span Time span where to solve the differential equation
 * @param sample_amount Desired amount of samples between the time span
 * @param initial_conditions Vector of initial conditions for the results
 * @param generators Vector of differential equations to be solved
 * @param point_solver Solver generator for the differential points
 * @return DifferentialEquationResult<T> Result for the multiple ODEs
 */
template<typename T = double>
DifferentialEquationResult<T> solve(const std::pair<double, double> &span,
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

  return DifferentialEquationResult(t, y);
}

/**
 * @brief Result from resolving an Nth order differential equation system
 * 
 * @tparam T Template data type for the differential equation system
 */
template<typename T = double>
struct NthOrderDifferentialEquationResult {
  /** Time vector for the results */
  std::vector<double> t;

  /**
   * @brief Vector of matrices to hold the results
   * Where:
   * Each vector in the vector represents an output result vector
   * and each point in the sub-vector represents a resolution point
   */
  std::vector<std::vector<T>> y;

  NthOrderDifferentialEquationResult(const std::vector<double> &_t, const std::vector<std::vector<T>> &_y) : t(_t), y(_y) {}

  NthOrderDifferentialEquationResult(std::vector<double> &&_t, std::vector<std::vector<T>> &&_y) : t(std::move(_t)), y(std::move(_y)) {}

  NthOrderDifferentialEquationResult(const NthOrderDifferentialEquationResult &other) : t(other.t), y(other.y) {}

  NthOrderDifferentialEquationResult(NthOrderDifferentialEquationResult &&other) : t(std::move(other.t)), y(std::move(other.y)) {}
};


template<typename T = double>
NthOrderDifferentialEquationResult<T> solveNthOrder(std::pair<double, double> span,
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
  std::vector<std::vector<T>> y;

  // Initial conditions for result matrix
  for (uint64_t eq_idx = 0; eq_idx < eq_amount; eq_idx++) {
    y.emplace_back(std::vector<double>(sample_amount));
    y.at(eq_idx).at(0) = initial_conditions.at(eq_idx);
  }

  // Now we need to solve the N differential equations of first order
  for (uint64_t t_idx = 1; t_idx < sample_amount; t_idx++) {
    t.at(t_idx) = span.first + (h * double(t_idx));

    for (uint64_t eq_idx = 0; eq_idx < eq_amount; eq_idx++) {
      y.at(eq_idx).at(t_idx) = point_solver(t.at(t_idx - 1), y.at(eq_idx).at(t_idx - 1), h, [&y, &A, &B, &t_idx, &eq_idx, &input_generator] (double t, const T &_y) {
        auto a_x = A.getRow(eq_idx) * _y;
        auto b_u = B.getRow(eq_idx) * input_generator(t);
        return a_x.at(0, 0) + b_u.at(0, 0);
      });
    }
  }

  return NthOrderDifferentialEquationResult(t, y);
}

};

#endif // __DIFFERENTIALECUATIONSOLVERS_HPP__