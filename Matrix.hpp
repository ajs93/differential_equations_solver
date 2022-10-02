#ifndef __MATRIX_HPP__
#define __MATRIX_HPP__

#include <string.h>
#include <stdint.h>

#include <sstream>
#include <stdexcept>
#include <vector>

/**
 * @brief Template class to define a generic matrix of N dimensions
 * 
 * @tparam T Template parameter type that will contain this matrix
 */
template<typename T = double>
class MatrixND {
public:
  /** Deleted dfault constructor */
  MatrixND() = delete;

  /**
   * @brief Matrix constructor
   * @exception std::runtime_error If any error is found in the dimensions vector
   * 
   * @param dimensions Vector of dimensions for this matrix
   */
  MatrixND(const std::vector<uint64_t> &dimensions) : _dimensions(dimensions), _values(_internalSize()) {
    for (const uint64_t &dim : _dimensions) {
      if (dim == 0) {
        throw std::runtime_error("Cannot create matrix with at least one dimension in zero");
      }
    }
  }

  /**
   * @brief Matrix constructor
   * @exception std::runtime_error If any error is found in the dimensions vector
   * 
   * @param dimensions Vector of dimensions for this matrix
   * @param initial_value Initial value for all of the indexes of this matrix
   */
  MatrixND(const std::vector<uint64_t> &dimensions, const T &initial_value) : _dimensions(dimensions), _values(_internalSize()) {
    for (const uint64_t &dim : _dimensions) {
      if (dim == 0) {
        throw std::runtime_error("Cannot create matrix with at least one dimension in zero");
      }
    }

    for (uint64_t idx = 0; idx < _values.size(); idx++) {
      _values.at(idx) = initial_value;
    }
  }

  /**
   * @brief Copy constructor from another N dimensional matrix
   * 
   * @param other Matrix from which to create the new one
   */
  MatrixND(const MatrixND &other) : _dimensions(other._dimensions), _values(other._values) {}

  /**
   * @brief Move constructor from another N dimensional matrix
   * 
   * @param other Matrix from which to create the new one
   */
  MatrixND(MatrixND &&other) : _dimensions(std::move(other._dimensions)), _values(std::move(other._values)) {}

  /** Default destructor */
  ~MatrixND() = default;

  /**
   * @brief Access a value of the matrix at @a coords
   * @exception std::runtime_error If any error is found in the coordinates to access the element
   * 
   * @param coords Coordinates inside the matrix of the element to be accesed
   * @return T& Element of the matrix at @a coords
   */
  T &at(const std::vector<uint64_t> &coords) {
    return _values.at(_indexOf(coords));
  }

  /**
   * @brief Access a value of the matrix at @a coords
   * @exception std::runtime_error If any error is found in the coordinates to access the element
   * 
   * @param coords Coordinates inside the matrix of the element to be accesed
   * @return T& Element of the matrix at @a coords
   */
  const T &at(const std::vector<uint64_t> &coords) const {
    return _values.at(_indexOf(coords));
  }

  /**
   * @brief Equality operator
   * 
   * @param other Matrix to be copied
   * @return MatrixND& Already copied matrix
   */
  MatrixND &operator=(const MatrixND &other) {
    _dimensions = other._dimensions;
    _values = other._values;
    return *this;
  }

  /**
   * @brief Add operator against another matrix
   * @exception std::runtime_error If the matrices cannot be added due to dimension issues
   * 
   * @param other Matrix to be added
   * @return MatrixND Matrix that results of the addition
   */
  MatrixND operator+(const MatrixND &other) const {
    if (!_sumSubDimensionCheck(other)) {
      std::ostringstream err;
      err << "Cannot sum matrices of dimensions: " << _printDimensions(_dimensions) << " and " << _printDimensions(other._dimensions);
      throw std::runtime_error(err.str());
    }

    MatrixND ret(*this);

    for (uint64_t idx = 0; idx < ret._values.size(); idx++) {
      ret._values.at(idx) += other._values.at(idx);
    }

    return ret;
  }

  /**
   * @brief Add and assign operator against another matrix
   * @exception std::runtime_error If the matrices cannot be added due to dimension issues
   * 
   * @param other Matrix to be added
   * @return MatrixND& Matrix that results of the addition
   */
  MatrixND &operator+=(const MatrixND &other) {
    if (!_sumSubDimensionCheck(other)) {
      std::ostringstream err;
      err << "Cannot sum matrices of dimensions: " << _printDimensions(_dimensions) << " and " << _printDimensions(other._dimensions);
      throw std::runtime_error(err.str());
    }

    for (uint64_t idx = 0; idx < _values.size(); idx++) {
      _values.at(idx) += other._values.at(idx);
    }

    return *this;
  }

  /**
   * @brief Add operator against a scalar
   * 
   * @param other Scalar to be added with the matrix
   * @return MatrixND Matrix that results of the addition
   */
  MatrixND operator+(const T &o) const {
    MatrixND ret(*this);

    for (uint64_t idx = 0; idx < ret._values.size(); idx++) {
      ret._values.at(idx) += o;
    }

    return ret;
  }

  /**
   * @brief Add operator against a scalar and an ND matrix
   * 
   * @param o Scalar to be added with the matrix
   * @param m Matrix to be added with the scalar
   * @return MatrixND Result of the addition
   */
  friend MatrixND operator+(const T &o, const MatrixND &m) {
    return m + o;
  }

  /**
   * @brief Add and assign operator against a scalar
   * 
   * @param o Scalar to be added with the matrix
   * @return MatrixND Matrix that results of the addition
   */
  MatrixND &operator+=(const T &o) {
    for (uint64_t idx = 0; idx < _values.size(); idx++) {
      _values.at(idx) += o;
    }

    return *this;
  }

  /**
   * @brief Substract operator against another matrix
   * @exception std::runtime_error If the matrices cannot be substracted due to dimension issues
   * 
   * @param other Matrix to be substracted
   * @return MatrixND Matrix that results of the substraction
   */
  MatrixND operator-(const MatrixND &other) const {
    if (!_sumSubDimensionCheck(other)) {
      std::ostringstream err;
      err << "Cannot substract matrices of dimensions: " << _printDimensions(_dimensions) << " and " << _printDimensions(other._dimensions);
      throw std::runtime_error(err.str());
    }

    MatrixND ret(*this);

    for (uint64_t idx = 0; idx < ret._values.size(); idx++) {
      ret._values.at(idx) -= other._values.at(idx);
    }

    return ret;
  }

  /**
   * @brief Substract and assign operator against another matrix
   * @exception std::runtime_error If the matrices cannot be substracted due to dimension issues
   * 
   * @param other Matrix to be substracted
   * @return MatrixND Matrix that results of the substraction
   */
  MatrixND &operator-=(const MatrixND &other) {
    if (!_sumSubDimensionCheck(other)) {
      std::ostringstream err;
      err << "Cannot substract matrices of dimensions: " << _printDimensions(_dimensions) << " and " << _printDimensions(other._dimensions);
      throw std::runtime_error(err.str());
    }

    for (uint64_t idx = 0; idx < _values.size(); idx++) {
      _values.at(idx) += other._values.at(idx);
    }

    return *this;
  }

  /**
   * @brief Substract operator against a scalar
   * 
   * @param other Scalar to be substracted
   * @return MatrixND Matrix that results of the substraction
   */
  MatrixND operator-(const T &o) const {
    MatrixND ret(*this);

    for (uint64_t idx = 0; idx < ret._values.size(); idx++) {
      ret._values.at(idx) -= o;
    }

    return ret;
  }

  /**
   * @brief Substract operator against a scalar and an ND matrix
   * 
   * @param o Scalar to be substracted with the matrix
   * @param m Matrix to be substracted with the scalar
   * @return MatrixND Result of the substraction
   */
  friend MatrixND operator-(const T &o, const MatrixND &m) {
    return m - o;
  }

  /**
   * @brief Substract and assign operator against a scalar
   * 
   * @param other Scalar to be substracted
   * @return MatrixND Matrix that results of the substraction
   */
  MatrixND &operator-=(const T &o) {
    for (uint64_t idx = 0; idx < _values.size(); idx++) {
      _values.at(idx) -= o;
    }

    return *this;
  }

  /**
   * @brief Multiply operator against a scalar
   * 
   * @param other Scalar to be multiplied
   * @return MatrixND Matrix that results of the multiplication
   */
  MatrixND operator*(const T &o) const {
    MatrixND ret(*this);

    for (uint64_t idx = 0; idx < ret._values.size(); idx++) {
      ret._values.at(idx) *= o;
    }

    return ret;
  }

  /**
   * @brief Multiply operator against a scalar and an ND matrix
   * 
   * @param o Scalar to be multiplicated with the matrix
   * @param m Matrix to be multiplicated with the scalar
   * @return MatrixND Result of the multiplication
   */
  friend MatrixND operator*(const T &o, const MatrixND &m) {
    return m * o;
  }

  /**
   * @brief Multiply and assign operator against a scalar
   * 
   * @param other Scalar to be multiplied
   * @return MatrixND Matrix that results of the multiplication
   */
  MatrixND &operator*=(const T &o) {
    for (uint64_t idx = 0; idx < _values.size(); idx++) {
      _values.at(idx) *= o;
    }

    return *this;
  }

  /**
   * @brief Divide operator against a scalar
   * 
   * @param other Scalar to be divided
   * @return MatrixND Matrix that results of the division
   */
  MatrixND operator/(const T &o) const {
    MatrixND ret(*this);

    for (uint64_t idx = 0; idx < ret._values.size(); idx++) {
      ret._values.at(idx) /= o;
    }

    return ret;
  }

  /**
   * @brief Divide operator against a scalar and an ND matrix
   * 
   * @param o Scalar to be divided with the matrix
   * @param m Matrix to be divided with the scalar
   * @return MatrixND Result of the division
   */
  friend MatrixND operator/(const T &o, const MatrixND &m) {
    return m / o;
  }

  /**
   * @brief Divide and assign operator against a scalar
   * 
   * @param other Scalar to be divided
   * @return MatrixND Matrix that results of the division
   */
  MatrixND &operator/=(const T &o) {
    for (uint64_t idx = 0; idx < _values.size(); idx++) {
      _values.at(idx) /= o;
    }

    return *this;
  }

  /**
   * @brief Get the dimensions of the matrix
   * 
   * @return const std::vector<uint64_t>& Vector containing the dimensions of this matrix
   */
  const std::vector<uint64_t> &getDimensions() const {
    return _dimensions;
  }

protected:
  /** Dimensions for this matrix */
  std::vector<uint64_t> _dimensions;

  /** Vector to hold the values of this matrix */
  std::vector<T> _values;

  /**
   * @brief Check if this matrix can be added with another matrix
   * 
   * @param other Matrix to be checked if can be added with this one
   * @return true Matrices can be added together
   * @return false Matrices can not be added together
   */
  bool _sumSubDimensionCheck(const MatrixND &other) const {
    for (uint64_t idx = 0; idx < _dimensions.size(); idx++) {
      if (_dimensions.at(idx) != other._dimensions.at(idx)) {
        return false;
      }
    }

    return true;
  }

  /**
   * @brief Get the internal size of the matrix
   * @note This denotes the amount of indexes in the matrix
   * 
   * @return uint64_t Internal size of the matrix
   */
  uint64_t _internalSize() const {
    uint64_t ret = 1;

    for (uint64_t idx = 0; idx < _dimensions.size(); idx++) {
      ret *= _dimensions[idx];
    }

    return ret;
  }

  /**
   * @brief Get the index in the value vector based on certain coordinates
   * @exception std::runtime_error If the coordinates are not valid
   * 
   * @param coords Coordinates to obtain the index
   * @return uint64_t Index in the value vector for the passed coordinates
   */
  uint64_t _indexOf(const std::vector<uint64_t> &coords) const {
    if (coords.size() != _dimensions.size()) {
      // NOTE: This shall never happen as we assure recourse creation is correct initialization
      // and this method is used in a protected context
      throw std::runtime_error("Invalid coordinates size");
    }

    for (uint64_t idx = 0; idx < coords.size(); idx++) {
      if (coords.at(idx) >= _dimensions.at(idx)) {
        std::ostringstream err;
        err << "Index " << _printDimensions(coords) << "} is out of range";
        throw std::runtime_error(err.str());
      }
    }

    uint64_t offset = 0;

    for (uint64_t idx = 0; idx < coords.size(); idx++) {
      uint64_t accum = 1;

      for (uint64_t accum_idx = idx + 1; accum_idx < coords.size(); accum_idx++) {
        accum *= _dimensions.at(accum_idx);
      }

      offset += coords.at(idx) * accum;
    }

    return offset;
  }

  /**
   * @brief Print the dimensions in a human-readable way
   * 
   * @param d Dimensions to be printed
   * @return std::string Dimensions in a human-readable way
   */
  static std::string _printDimensions(const std::vector<uint64_t> &d) {
    std::ostringstream out;
    out << "{";

    for (uint64_t idx = 0; idx < d.size(); idx++) {
      out << d.at(idx);

      if (idx != d.size() - 1) {
        out << "x";
      }
    }

    out << "}";

    return out.str();
  }
};

/**
 * @brief Two dimensional matrix specialization
 * 
 * @tparam T Template parameter type that will contain this matrix
 */
template<typename T = double>
class Matrix2D : public MatrixND<T> {
public:
  using MatrixND<T>::_indexOf;
  using MatrixND<T>::_printDimensions;

  /** Deleted default constructor */
  Matrix2D() = delete;

  /**
   * @brief 2D Matrix constructor
   * 
   * @param rows Amount of rows for the matrix
   * @param columns Amount of columns for the matrix
   */
  Matrix2D(uint64_t rows, uint64_t columns) : MatrixND<T>({rows, columns}), __rows(rows), __columns(columns) {}

  /**
   * @brief 2D Matrix constructor with an initial value
   * 
   * @param rows Amount of rows for the matrix
   * @param columns Amount of columns for the matrix
   * @param initial_value Initial value for all indexes in the matrix
   */
  Matrix2D(uint64_t rows, uint64_t columns, const T &initial_value) : MatrixND<T>({rows, columns}, initial_value), __rows(rows), __columns(columns) {}

  /**
   * @brief 2D Matrix copy constructor
   * 
   * @param other Matrix to be copied from
   */
  Matrix2D(const Matrix2D &other) : MatrixND<T>(other), __rows(other.__rows), __columns(other.__columns) {}

  /**
   * @brief 2D Matrix move constructor
   * 
   * @param other Matrix to be moved from
   */
  Matrix2D(Matrix2D &&other) : MatrixND<T>(std::move(other)), __rows(other.__rows), __columns(other.__columns) {}

  /**
   * @brief 2D Matrix constructor from row vectors
   * @note This constructor is very usefull to create matrices quickly in code
   * 
   * @exception std::runtime_error If the rows are not consistent
   * 
   * @param rows Vector of vectors that represent the different rows of the desired 2D matrix
   * @note The row and column amount will be deduced from the @a rows argument
   */
  Matrix2D(const std::vector<std::vector<T>> &rows) : MatrixND<T>({rows.size(), rows.at(0).size()}), __rows(rows.size()), __columns(rows.at(0).size()) {
    for (const auto &row : rows) {
      if (row.size() != __columns) {
        throw std::runtime_error("Cannot create two dimensional matrix if rows sizes doesn't match");
      }
    }

    for (uint64_t r = 0; r < __rows; r++) {
      for (uint64_t c = 0; c < __columns; c++) {
        this->at(r, c) = rows.at(r).at(c);
      }
    }
  }

  /** Default destructor */
  ~Matrix2D() = default;

  /**
   * @brief Equality operator
   * 
   * @param other Matrix to be copied inside this object
   * @return Matrix2D& Same matrix reference
   */
  Matrix2D &operator=(const Matrix2D &other) {
    this->_values = other._values;
    this->_dimensions = other._dimensions;
    this->__rows = other.__rows;
    this->__columns = other.__columns;
    return *this;
  }

  /**
   * @brief Access to an item in the matrix
   * 
   * @param row Row index
   * @param column Column index
   * @return T& Desired item
   */
  T &at(uint64_t row, uint64_t column) {
    return this->_values.at(_indexOf({row, column}));
  }

  /**
   * @brief Access to an item in the matrix
   * 
   * @param row Row index
   * @param column Column index
   * @return T& Desired item
   */
  const T &at(uint64_t row, uint64_t column) const {
    return this->_values.at(_indexOf({row, column}));
  }

  /**
   * @brief Add operator against another matrix
   * @exception std::runtime_error If the matrices cannot be added due to dimension issues
   * 
   * @param other Matrix to be added
   * @return Matrix2D Matrix that results of the addition
   */
  Matrix2D operator+(const Matrix2D &other) const {
    auto ret = MatrixND<T>::operator+(other);
    return Matrix2D::fromMatrixND(ret);
  }

  /**
   * @brief Add operator against a scalar
   * 
   * @param other Scalar to be added with the matrix
   * @return Matrix2D Matrix that results of the addition
   */
  Matrix2D operator+(const T &o) const {
    Matrix2D ret(*this);
    ret += o;
    return ret;
  }

  /**
   * @brief Addition operator against a scalar and a 2D matrix
   * 
   * @param o Scalar to be added with the matrix
   * @param m Matrix to be added with the scalar
   * @return MatrixND Result of the addition
   */
  friend Matrix2D operator+(const T &o, const Matrix2D &m) {
    return m + o;
  }

  /**
   * @brief Substract operator against another matrix
   * @exception std::runtime_error If the matrices cannot be substracted due to dimension issues
   * 
   * @param other Matrix to be substracted
   * @return Matrix2D Matrix that results of the Substraction
   */
  Matrix2D operator-(const Matrix2D &other) const {
    auto ret = MatrixND<T>::operator-(other);
    return Matrix2D::fromMatrixND(ret);
  }

  /**
   * @brief Substract operator against a scalar
   * 
   * @param other Scalar to be substracted with the matrix
   * @return Matrix2D Matrix that results of the substraction
   */
  Matrix2D operator-(const T &o) const {
    Matrix2D ret(*this);
    ret -= o;
    return ret;
  }

  /**
   * @brief Substract operator against a scalar and a 2D matrix
   * 
   * @param o Scalar to be substracted with the matrix
   * @param m Matrix to be substracted with the scalar
   * @return MatrixND Result of the substraction
   */
  friend Matrix2D operator-(const T &o, const Matrix2D &m) {
    return m - o;
  }

  /**
   * @brief Multiply operator against a scalar
   * 
   * @param other Scalar to be multiplied with the matrix
   * @return Matrix2D Matrix that results of the multiplication
   */
  Matrix2D operator*(const T &o) const {
    Matrix2D ret(*this);
    ret *= o;
    return ret;
  }

  /**
   * @brief Multiply operator against a scalar and a 2D matrix
   * 
   * @param o Scalar to be multiplied with the matrix
   * @param m Matrix to be multiplied with the scalar
   * @return MatrixND Result of the multiplication
   */
  friend Matrix2D operator*(const T &o, const Matrix2D &m) {
    return m * o;
  }

  /**
   * @brief Divide operator against a scalar
   * 
   * @param other Scalar to be divided with the matrix
   * @return Matrix2D Matrix that results of the division
   */
  Matrix2D operator/(const T &o) const {
    Matrix2D ret(*this);
    ret /= o;
    return ret;
  }

  /**
   * @brief Divide operator against a scalar and a 2D matrix
   * 
   * @param o Scalar to be divided with the matrix
   * @param m Matrix to be divided with the scalar
   * @return MatrixND Result of the division
   */
  friend Matrix2D operator/(const T &o, const Matrix2D &m) {
    return m / o;
  }

  /**
   * @brief Multiply operator with another matrix
   * @exception std::runtime_error If the two matrices cannot be multiplied due to their respective dimensions
   * 
   * @param other Matrix to be multiplied with
   * @return Matrix2D Matrix that results of the multiplication
   */
  Matrix2D operator*(const Matrix2D &other) const {
    if (__columns != other.__rows) {
      std::ostringstream err;
      err << "Cannot multiply 2D matrixes of " << _printDimensions(this->_dimensions) << " and " << other._printDimensions(other._dimensions);
      throw std::runtime_error(err.str());
    }

    Matrix2D ret(__rows, other.__columns);

    for (uint64_t ridx = 0; ridx < ret.__rows; ridx++) {
      for (uint64_t cidx = 0; cidx < ret.__columns; cidx++) {
        T v(0);

        for (uint64_t nidx = 0; nidx < other.__rows; nidx++) {
          v += this->at(ridx, nidx) * other.at(nidx, cidx);
        }

        ret.at(ridx, cidx) = v;
      }
    }

    return ret;
  }

  /**
   * @brief Get a row from the matrix
   * @exception std::runtime_error If the row index is out of range
   * 
   * @param row_idx Desired row index
   * @return Matrix2D Desired row
   */
  Matrix2D getRow(uint64_t row_idx) const {
    if (row_idx >= this->__rows) {
      std::ostringstream err;
      err << "Row index " << row_idx << " is out of range (" << this->__rows << ")";
      throw std::runtime_error(err.str());
    }
    
    Matrix2D ret(1, this->__columns);

    for (uint64_t idx = 0; idx < this->__columns; idx++) {
      ret.at(0, idx) = this->at(row_idx, idx);
    }

    return ret;
  }

  /**
   * @brief Get a column from the matrix
   * @exception std::runtime_error If the column index is out of range
   * 
   * @param row_idx Desired column index
   * @return Matrix2D Desired column
   */
  Matrix2D getColumn(uint64_t col_idx) const {
    if (col_idx >= this->__columns) {
      std::ostringstream err;
      err << "Column index " << col_idx << " is out of range (" << this->__columns << ")";
      throw std::runtime_error(err.str());
    }
    
    Matrix2D ret(this->__rows, 1);

    for (uint64_t idx = 0; idx < this->__rows; idx++) {
      ret.at(idx, 0) = this->at(idx, col_idx);
    }

    return ret;
  }

  /**
   * @brief Get the row amount of this matrix
   * 
   * @return uint64_t Row amount
   */
  uint64_t getRowAmount() const {
    return __rows;
  }

  /**
   * @brief Get the column amount of this matrix
   * 
   * @return uint64_t Column amount
   */
  uint64_t getColumnAmount() const {
    return __columns;
  }

  /**
   * @brief Get the total size of this matrix
   * 
   * @return std::size_t Total size of the matrix
   */
  std::size_t getSize() const {
    return this->_values.size();
  }

  /**
   * @brief Construct a 2D Matrix from an N dimensional one
   * @exception std::runtime_error If the N dimensional matrix is not suited to create
   * a 2D one
   * 
   * @param other Matrix to be used as base
   * @return Matrix2D Resulting matrix from the transformation
   */
  static Matrix2D fromMatrixND(const MatrixND<T> &other) {
    if (other.getDimensions().size() != 2) {
      std::ostringstream err;
      err << "Cannot create 2D matrix from " << other.getDimensions().size() << " dimensional matrix";
      throw std::runtime_error(err.str());
    }

    Matrix2D ret(other.getDimensions().at(0), other.getDimensions().at(1));
    ret._values = static_cast<const Matrix2D*>(&other)->_values;
    return ret;
  }

private:
  /** Row amount of this matrix */
  uint64_t __rows;

  /** Column amount of this matrix */
  uint64_t __columns;
};

#endif /* __MATRIX_HPP__ */