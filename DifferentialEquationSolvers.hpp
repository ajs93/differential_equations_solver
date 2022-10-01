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
 * @brief Result from resolving an Nth order differential equation system
 * 
 * @tparam T Template data type for the differential equation system
 */
template<typename T = double>
struct DifferentialEquationResult {
  /** Time vector for the results */
  std::vector<double> t;

  /**
   * @brief Vector of matrices to hold the results
   * Where:
   * Each vector in the vector represents an output result vector
   * and each point in the sub-vector represents a resolution point
   */
  std::vector<std::vector<T>> y;

  DifferentialEquationResult(const std::vector<double> &_t, const std::vector<std::vector<T>> &_y) : t(_t), y(_y) {}

  DifferentialEquationResult(std::vector<double> &&_t, std::vector<std::vector<T>> &&_y) : t(std::move(_t)), y(std::move(_y)) {}

  DifferentialEquationResult(const DifferentialEquationResult &other) : t(other.t), y(other.y) {}

  DifferentialEquationResult(DifferentialEquationResult &&other) : t(std::move(other.t)), y(std::move(other.y)) {}

  DifferentialEquationResult &operator=(const DifferentialEquationResult &other) {
    t = other.t;
    y = other.y;
    return *this;
  }
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
  std::vector<std::vector<T>> y;

  t.at(0) = span.first;

  for (uint64_t gen_idx = 0; gen_idx < generators.size(); gen_idx++) {
    y.emplace_back(std::vector<double>(sample_amount));
    y.at(gen_idx).at(0) = initial_conditions.at(gen_idx);
  }

  for (uint64_t t_idx = 1; t_idx < sample_amount; t_idx++) {
    t.at(t_idx) = span.first + (h * double(t_idx));

    for (uint64_t gen_idx = 0; gen_idx < generators.size(); gen_idx++) {
      y.at(gen_idx).at(t_idx) = point_solver(t.at(t_idx - 1), y.at(gen_idx).at(t_idx - 1), h, generators.at(gen_idx));
    }
  }

  return DifferentialEquationResult(t, y);
}


template<typename T = double>
DifferentialEquationResult<T> solveNthOrder(std::pair<double, double> span,
                                            uint64_t sample_amount,
                                            const Matrix2D<T> &initial_conditions,
                                            const InputGenerator<T> &input_generator,
                                            const DiffPointSolver<T> &point_solver,
                                            const Matrix2D<T> &A,
                                            const Matrix2D<T> &B,
                                            const Matrix2D<T> &C,
                                            const Matrix2D<T> &D) {
  if (A.getRowAmount() != A.getColumnAmount()) {
    throw std::runtime_error("A matrix must be a square matrix");
  }

  if (B.getColumnAmount() != 1) {
    throw std::runtime_error("B matrix must be a column matrix");
  }

  if (A.getRowAmount() != B.getRowAmount()) {
    throw std::runtime_error("A and B matrices must have the same row amount");
  }

  if (C.getColumnAmount() != A.getColumnAmount()) {
    throw std::runtime_error("A and C matrices must have the same column amount");
  }

  if (C.getRowAmount() != D.getRowAmount()) {
    throw std::runtime_error("C and D matrices must have the same row amount");
  }

  const double h = (span.second - span.first) / double(sample_amount - 1);
  std::vector<double> t(sample_amount);
  std::vector<std::vector<T>> y;

  Matrix2D<T> xp = initial_conditions; // Memory allocation optimization
  Matrix2D<T> aux_xp(xp.getRowAmount(), xp.getColumnAmount()); // Memory allocation optimization

  for (uint64_t idx = 0; idx < C.getRowAmount(); idx++) {
    y.emplace_back(std::vector<T>(sample_amount));
    y.at(idx).at(0) = ((C.getRow(idx) * initial_conditions) + (D.getRow(idx) * input_generator(0))).at(0, 0);
  }

  for (uint64_t t_idx = 1; t_idx < sample_amount; t_idx++) {
    t.at(t_idx) = span.first + (h * double(t_idx));

    for (uint64_t xp_row = 0; xp_row < xp.getRowAmount(); xp_row++) {
      aux_xp.at(xp_row, 0) = point_solver(t.at(t_idx - 1), xp.at(xp_row, 0), h, [&A, &xp, &B, &input_generator, &xp_row, &t_idx](double t, const T &_y) {
        return (A.getRow(xp_row) * xp).at(0, 0) + B.at(xp_row, 0) * input_generator(t);
      });
    }
    
    for (uint64_t idx = 0; idx < y.size(); idx++) {
      y.at(idx).at(t_idx) = ((C.getRow(idx) * aux_xp) + (D.getRow(idx) * input_generator(t.at(t_idx)))).at(0, 0);
    }

    xp = aux_xp;
  }

  return DifferentialEquationResult(t, y);
}

};

#endif // __DIFFERENTIALECUATIONSOLVERS_HPP__