#include <catch2/catch_test_macros.hpp>

#include <catch2/catch_approx.hpp>

#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>

#include <Matrix.hpp>

SCENARIO("Two dimensional matrix constructors", "[auto]") {
  GIVEN("Several valid rows and columns amount and several initializer values") {
    uint64_t column_amount = GENERATE(take(10, random(1, 20)));
    uint64_t row_amount = GENERATE(take(10, random(1, 20)));
    double initializer_value = GENERATE(take(10, random(-100.0, 100.0)));

    WHEN("A matrix is created with row amount in zero and any other column amount") {
      THEN("A runtime exception shall be thrown") {
        REQUIRE_THROWS_AS(Matrix2D(0, column_amount), std::runtime_error);
        REQUIRE_THROWS_AS(Matrix2D(0, column_amount, initializer_value), std::runtime_error);
      }
    }

    WHEN("A matrix is created with any row amount and the column amount in zero") {
      THEN("A runtime exception shall be thrown") {
        REQUIRE_THROWS_AS(Matrix2D(row_amount, 0), std::runtime_error);
        REQUIRE_THROWS_AS(Matrix2D(row_amount, 0, initializer_value), std::runtime_error);
      }
    }
  }

  AND_GIVEN("Two valid row initializers with two columns with known values") {
    std::vector<std::vector<double>> rows = {
      {0, 1},
      {2, 3},
    };

    WHEN("A matrix is created with the row initializers") {
      Matrix2D m(rows);

      THEN("The row amount is correct") {
        REQUIRE(m.getRowAmount() == 2);
      }

      AND_THEN("The column amount is correct") {
        REQUIRE(m.getColumnAmount() == 2);
      }

      AND_THEN("The matrix size is correct") {
        REQUIRE(m.getSize() == 4);
      }

      AND_THEN("The indexes are correct") {
        REQUIRE(m.at(0, 0) == 0);
        REQUIRE(m.at(0, 1) == 1);
        REQUIRE(m.at(1, 0) == 2);
        REQUIRE(m.at(1, 1) == 3);
      }
    }
  }

  AND_GIVEN("Two valid row initializers with four columns with known values") {
    std::vector<std::vector<double>> rows = {
      {0, 1, 2, 3},
      {4, 5, 6, 7},
    };

    WHEN("A matrix is created with the row initializers") {
      Matrix2D m(rows);

      THEN("The row amount is correct") {
        REQUIRE(m.getRowAmount() == 2);
      }

      AND_THEN("The column amount is correct") {
        REQUIRE(m.getColumnAmount() == 4);
      }

      AND_THEN("The matrix size is correct") {
        REQUIRE(m.getSize() == 8);
      }

      AND_THEN("The indexes are correct") {
        REQUIRE(m.at(0, 0) == 0);
        REQUIRE(m.at(0, 1) == 1);
        REQUIRE(m.at(0, 2) == 2);
        REQUIRE(m.at(0, 3) == 3);
        REQUIRE(m.at(1, 0) == 4);
        REQUIRE(m.at(1, 1) == 5);
        REQUIRE(m.at(1, 2) == 6);
        REQUIRE(m.at(1, 3) == 7);
      }
    }
  }

  AND_GIVEN("Four valid row initializers with two columns with known values") {
    std::vector<std::vector<double>> rows = {
      {0, 1},
      {2, 3},
      {4, 5},
      {6, 7},
    };

    WHEN("A matrix is created with the row initializers") {
      Matrix2D m(rows);

      THEN("The row amount is correct") {
        REQUIRE(m.getRowAmount() == 4);
      }

      AND_THEN("The column amount is correct") {
        REQUIRE(m.getColumnAmount() == 2);
      }

      AND_THEN("The matrix size is correct") {
        REQUIRE(m.getSize() == 8);
      }

      AND_THEN("The indexes are correct") {
        REQUIRE(m.at(0, 0) == 0);
        REQUIRE(m.at(0, 1) == 1);
        REQUIRE(m.at(1, 0) == 2);
        REQUIRE(m.at(1, 1) == 3);
        REQUIRE(m.at(2, 0) == 4);
        REQUIRE(m.at(2, 1) == 5);
        REQUIRE(m.at(3, 0) == 6);
        REQUIRE(m.at(3, 1) == 7);
      }
    }
  }

  AND_GIVEN("Two invalid row initializers with different column sizes") {
    std::vector<std::vector<double>> rows = {
      {0, 1},
      {2, 3, 4},
    };

    WHEN("A matrix is created with the row initializers") {
      THEN("A runtime exception shall be thrown") {
        REQUIRE_THROWS_AS(Matrix2D(rows), std::runtime_error);
      }
    }
  }

  AND_GIVEN("An N-dimensional matrix of two dimensions") {
    uint64_t column_amount = GENERATE(take(10, random(1, 20)));
    uint64_t row_amount = GENERATE(take(10, random(1, 20)));
    MatrixND ndim({row_amount, column_amount});

    WHEN("A matrix is created from the N-dimensional matrix") {
      Matrix2D m = Matrix2D<double>::fromMatrixND(ndim);

      THEN("The created 2D matrix shall have valid row and column amount") {
        REQUIRE(m.getRowAmount() == row_amount);
        REQUIRE(m.getColumnAmount() == column_amount);
      }
    }
  }

  AND_GIVEN("An N-dimensional matrix of more than two dimensions (all different from zero)") {
    uint64_t n_dimensions = GENERATE(take(10, random(3, 10)));
    std::vector<uint64_t> dimensions(n_dimensions);

    for (uint64_t idx = 0; idx < dimensions.size(); idx++) {
      dimensions.at(idx) = 3;
    }

    MatrixND ndim(dimensions);

    WHEN("A matrix is created from the N-dimensional matrix") {
      THEN("A runtime exception shall be thrown") {
        REQUIRE_THROWS_AS(Matrix2D<double>::fromMatrixND(ndim), std::runtime_error);
      }
    }
  }
}

SCENARIO("Matrix access methods", "[auto]") {
  GIVEN("A matrix of 1x1 with no initializer value") {
    Matrix2D m(1, 1);

    WHEN("Matrix value is accessed") {
      double val = m.at(0, 0);
      THEN("The value is zero") {
        REQUIRE(val == 0.0);
      }
    }
  }

  GIVEN("Matrix with valid rows and columns amount") {
    uint64_t row_amount = GENERATE(take(10, random(1, 10)));
    uint64_t column_amount = GENERATE(take(10, random(1, 10)));

    Matrix2D m(row_amount, column_amount);

    WHEN("Matrix characteristics are accessed") {
      THEN("The row amount is correct") {
        REQUIRE(m.getRowAmount() == row_amount);
      }

      AND_THEN("The column amount is correct") {
        REQUIRE(m.getColumnAmount() == column_amount);
      }

      AND_THEN("The matrix size is correct") {
        REQUIRE(m.getSize() == (row_amount * column_amount));
      }
    }

    WHEN("Valid indexes are accessed") {
      THEN("The indexes are initialized in zero") {
        for (uint64_t r = 0; r < row_amount; r++) {
          for (uint64_t c = 0; c < column_amount; c++) {
            REQUIRE(m.at(r, c) == 0.0);
          }
        }
      }
    }

    WHEN("Invalid indexes are accessed") {
      THEN("A runtime exception shall be thrown") {
        REQUIRE_THROWS_AS(m.at(row_amount + 1, 0), std::runtime_error);
        REQUIRE_THROWS_AS(m.at(0, column_amount + 1), std::runtime_error);
      }
    }

    WHEN("Valid rows are retrieved") {
      std::vector<Matrix2D<double>> rows;

      for (uint64_t r = 0; r < row_amount; r++) {
        rows.emplace_back(m.getRow(r));
      }

      THEN("The indexes are initialized in zero") {
        for (const auto &r : rows) {
          for (uint64_t c = 0; c < column_amount; c++) {
            REQUIRE(r.at(0, c) == 0.0);
          }
        }
      }
    }

    WHEN("Valid columns are retrieved") {
      std::vector<Matrix2D<double>> cols;

      for (uint64_t c = 0; c < column_amount; c++) {
        cols.emplace_back(m.getColumn(c));
      }

      THEN("The indexes are initialized in zero") {
        for (const auto &c : cols) {
          for (uint64_t r = 0; r < row_amount; r++) {
            REQUIRE(c.at(r, 0) == 0.0);
          }
        }
      }
    }

    WHEN("Invalid rows are retrieved") {
      THEN("A runtime exception shall be thrown") {
        REQUIRE_THROWS_AS(m.getRow(row_amount + 1), std::runtime_error);
      }
    }
    
    WHEN("Invalid columns are retrieved") {
      THEN("A runtime exception shall be thrown") {
        REQUIRE_THROWS_AS(m.getColumn(column_amount + 1), std::runtime_error);
      }
    }
  }

  GIVEN("Matrix with valid rows and columns amount and an initializer value") {
    uint64_t row_amount = GENERATE(take(10, random(1, 10)));
    uint64_t column_amount = GENERATE(take(10, random(1, 10)));
    double initializer_value = GENERATE(take(10, random(-100.0, 100.0)));

    Matrix2D m(row_amount, column_amount, initializer_value);

    WHEN("Matrix characteristics are accessed") {
      THEN("The row amount is correct") {
        REQUIRE(m.getRowAmount() == row_amount);
      }

      AND_THEN("The column amount is correct") {
        REQUIRE(m.getColumnAmount() == column_amount);
      }

      AND_THEN("The matrix size is correct") {
        REQUIRE(m.getSize() == (row_amount * column_amount));
      }
    }

    WHEN("Valid indexes are accessed") {
      THEN("The indexes are initialized in the value of the initializer value") {
        for (uint64_t r = 0; r < row_amount; r++) {
          for (uint64_t c = 0; c < column_amount; c++) {
            REQUIRE(m.at(r, c) == Catch::Approx(initializer_value).epsilon(0.01));
          }
        }
      }
    }

    WHEN("Invalid indexes are accessed") {
      THEN("A runtime exception shall be thrown") {
        REQUIRE_THROWS_AS(m.at(row_amount + 1, 0), std::runtime_error);
        REQUIRE_THROWS_AS(m.at(0, column_amount + 1), std::runtime_error);
      }
    }
  }
}

SCENARIO("Matrix operators", "[auto]") {
  GIVEN("A matrix with no initializer value") {
    uint64_t row_amount = GENERATE(take(10, random(1, 10)));
    uint64_t column_amount = GENERATE(take(10, random(1, 10)));
    Matrix2D m0(row_amount, column_amount);

    WHEN("The matrix is added a simple scalar") {
      double adder = GENERATE(take(10, random(-100.0, 100.0)));

      Matrix2D m1 = m0 + adder;
      
      THEN("The matrix shall be updated in all of it's indexes") {
        for (uint64_t r = 0; r < row_amount; r++) {
          for (uint64_t c = 0; c < column_amount; c++) {
            REQUIRE(m1.at(r, c) == Catch::Approx(adder).epsilon(0.01));
          }
        }
      }
    }

    AND_WHEN("The matrix is substracted a simple scalar") {
      double subs = GENERATE(take(10, random(-100.0, 100.0)));

      Matrix2D m1 = m0 - subs;
      
      THEN("The matrix shall be updated in all of it's indexes") {
        for (uint64_t r = 0; r < row_amount; r++) {
          for (uint64_t c = 0; c < column_amount; c++) {
            REQUIRE(m1.at(r, c) == Catch::Approx(subs * -1).epsilon(0.01));
          }
        }
      }
    }

    AND_WHEN("The matrix is added by one and multiplied by a simple scalar") {
      double mul = GENERATE(take(10, random(-100.0, 100.0)));

      Matrix2D m1 = (m0 + 1.0) * mul;
      
      THEN("The matrix shall be updated in all of it's indexes") {
        for (uint64_t r = 0; r < row_amount; r++) {
          for (uint64_t c = 0; c < column_amount; c++) {
            REQUIRE(m1.at(r, c) == Catch::Approx(mul).epsilon(0.01));
          }
        }
      }
    }

    AND_WHEN("The matrix is added by one and divided by a simple scalar different than zero") {
      double div = GENERATE(take(10, filter([](double d) { return d != 0.0; }, random(-100.0, 100.0))));

      Matrix2D m1 = (m0 + 1) / div;
      
      THEN("The matrix shall be updated in all of it's indexes") {
        for (uint64_t r = 0; r < row_amount; r++) {
          for (uint64_t c = 0; c < column_amount; c++) {
            REQUIRE(m1.at(r, c) == Catch::Approx(1 / div).epsilon(0.01));
          }
        }
      }
    }

    AND_WHEN("The matrix is added by another matrix of the same size with an initializer value") {
      double initializer_value = GENERATE(take(10, random(-100.0, 100.0)));

      Matrix2D other_matrix(row_amount, column_amount, initializer_value);
      Matrix2D m1 = m0 + other_matrix;

      THEN("The matrix shall be updated in all of it's indexes") {
        for (uint64_t r = 0; r < row_amount; r++) {
          for (uint64_t c = 0; c < column_amount; c++) {
            REQUIRE(m1.at(r, c) == Catch::Approx(other_matrix.at(r, c)).epsilon(0.01));
          }
        }
      }
    }

    AND_WHEN("The matrix is substracted by another matrix of the same size with an initializer value") {
      double initializer_value = GENERATE(take(10, random(-100.0, 100.0)));

      Matrix2D other_matrix(row_amount, column_amount, initializer_value);
      Matrix2D m1 = m0 - other_matrix;

      THEN("The matrix shall be updated in all of it's indexes") {
        for (uint64_t r = 0; r < row_amount; r++) {
          for (uint64_t c = 0; c < column_amount; c++) {
            REQUIRE(m1.at(r, c) == Catch::Approx(-1 * other_matrix.at(r, c)).epsilon(0.01));
          }
        }
      }
    }

    AND_WHEN("The matrix is added by another matrix of different size") {
      Matrix2D other_matrix_1(row_amount, column_amount + 1);
      Matrix2D other_matrix_2(row_amount + 1, column_amount);
      
      THEN("An exception shall be thrown") {
        REQUIRE_THROWS(m0 + other_matrix_1);
        REQUIRE_THROWS(m0 + other_matrix_2);
      }
    }

    AND_WHEN("The matrix is substracted by another matrix of different size") {
      Matrix2D other_matrix_1(row_amount, column_amount + 1);
      Matrix2D other_matrix_2(row_amount + 1, column_amount);
      
      THEN("An exception shall be thrown") {
        REQUIRE_THROWS(m0 - other_matrix_1);
        REQUIRE_THROWS(m0 - other_matrix_2);
      }
    }

    AND_WHEN("The matrix (with sizes MxN) is multiplied by another matrix with an initializer value (with sizes NxP)") {
      uint64_t p_dim = GENERATE(take(10, random(1, 10)));
      double initializer_value = GENERATE(take(5, random(-100.0, 100.0)));

      Matrix2D other_matrix(column_amount, p_dim, initializer_value);

      Matrix2D m1 = m0 * other_matrix;

      THEN("The resulting matrix shall be of dimensions MxP") {
        REQUIRE(m1.getRowAmount() == row_amount);
        REQUIRE(m1.getColumnAmount() == p_dim);
      }

      AND_THEN("The resulting matrix shall contain all of it's values in zero") {
        for (uint64_t r = 0; r < row_amount; r++) {
          for (uint64_t c = 0; c < p_dim; c++) {
            REQUIRE(m1.at(r, c) == 0.0);
          }
        }
      }
    }

    AND_WHEN("The matrix (with sizes MxN) is multiplied by another matrix (with sizes PxQ) with N != P") {
      static uint64_t ca;

      if (ca != column_amount) {
        ca = column_amount;
      }

      uint64_t p_dim = GENERATE(take(10, filter([=] (uint64_t dim) { return dim != ca; }, random(1, 10))));
      uint64_t q_dim = GENERATE(take(10, random(1, 10)));

      Matrix2D other_matrix(p_dim, q_dim);

      INFO("N dim: " << column_amount << " P dim: " << p_dim);

      THEN("An exception shall be thrown") {
        REQUIRE_THROWS(m0 * other_matrix);
      }
    }
  }
}

SCENARIO("Basic matrix multiplication checks", "[auto]") {
  GIVEN("Two 2x2 matrices with known initialized indexes") {
    Matrix2D m0({{45, 22}, {0, -8}});
    Matrix2D m1({{13, 8}, {-9, 5}});

    WHEN("The two matrices are multiplied") {
      Matrix2D res = m0 * m1;

      THEN("The dimensions of the result is 2x2") {
        REQUIRE(res.getRowAmount() == 2);
        REQUIRE(res.getColumnAmount() == 2);
      }

      AND_THEN("The result is well known") {
        REQUIRE(res.at(0, 0) == Catch::Approx(387).epsilon(0.01));
        REQUIRE(res.at(0, 1) == Catch::Approx(470).epsilon(0.01));
        REQUIRE(res.at(1, 0) == Catch::Approx(72).epsilon(0.01));
        REQUIRE(res.at(1, 1) == Catch::Approx(-40).epsilon(0.01));
      }
    }
  }

  GIVEN("One 3x2 matrix and another 2x3 matrix with known initialized indexes") {
    Matrix2D m0({{45, 22}, {0, -8}, {5, 6}});
    Matrix2D m1({{13, 8, 1}, {-9, 5, 0}});

    WHEN("The two matrices are multiplied") {
      Matrix2D res = m0 * m1;

      THEN("The dimensions of the result is 3x3") {
        REQUIRE(res.getRowAmount() == 3);
        REQUIRE(res.getColumnAmount() == 3);
      }

      AND_THEN("The result is well known") {
        REQUIRE(res.at(0, 0) == Catch::Approx(387).epsilon(0.01));
        REQUIRE(res.at(0, 1) == Catch::Approx(470).epsilon(0.01));
        REQUIRE(res.at(0, 2) == Catch::Approx(45).epsilon(0.01));
        REQUIRE(res.at(1, 0) == Catch::Approx(72).epsilon(0.01));
        REQUIRE(res.at(1, 1) == Catch::Approx(-40).epsilon(0.01));
        REQUIRE(res.at(1, 2) == Catch::Approx(0).epsilon(0.01));
        REQUIRE(res.at(2, 0) == Catch::Approx(11).epsilon(0.01));
        REQUIRE(res.at(2, 1) == Catch::Approx(70).epsilon(0.01));
        REQUIRE(res.at(2, 2) == Catch::Approx(5).epsilon(0.01));
      }
    }
  }

  GIVEN("One 1x3 matrix and another 3x1 matrix with known initialized indexes") {
    Matrix2D m0({{6, 11, -5}});
    Matrix2D m1(std::vector<std::vector<double>>{{0}, {5}, {-7}});

    WHEN("The two matrices are multiplied") {
      Matrix2D res = m0 * m1;

      THEN("The dimensions of the result is 1x1") {
        REQUIRE(res.getRowAmount() == 1);
        REQUIRE(res.getColumnAmount() == 1);
      }

      AND_THEN("The result is well known") {
        REQUIRE(res.at(0, 0) == Catch::Approx(90).epsilon(0.01));
      }
    }
  }
}