#ifndef __MATRIX_HPP__
#define __MATRIX_HPP__

#include <string.h>
#include <stdint.h>

#include <sstream>
#include <stdexcept>
#include <vector>

template<typename T = double>
class MatrixND {
public:
  MatrixND() = delete;

  MatrixND(const std::vector<uint64_t> &dimensions) : _dimensions(dimensions), _values(_internalSize()) {
    for (const uint64_t &dim : _dimensions) {
      if (dim == 0) {
        throw std::runtime_error("Cannot create matrix with at least one dimension in zero");
      }
    }
  }

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

  MatrixND(const MatrixND &other) : _dimensions(other._dimensions), _values(other._values) {}

  MatrixND(MatrixND &&other) : _dimensions(std::move(other._dimensions)), _values(std::move(other._values)) {}

  ~MatrixND() = default;

  T &at(const std::vector<uint64_t> &coords) {
    return _values.at(_indexOf(coords));
  }

  const T &at(const std::vector<uint64_t> &coords) const {
    return _values.at(_indexOf(coords));
  }

  MatrixND &operator=(const MatrixND &other) {
    _dimensions = other._dimensions;
    _values = other._values;
    return *this;
  }

  MatrixND operator+(const MatrixND &other) const {
    if (_sumSubDimensionCheck(other)) {
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

  MatrixND &operator+=(const MatrixND &other) {
    if (_sumSubDimensionCheck(other)) {
      std::ostringstream err;
      err << "Cannot sum matrices of dimensions: " << _printDimensions(_dimensions) << " and " << _printDimensions(other._dimensions);
      throw std::runtime_error(err.str());
    }

    for (uint64_t idx = 0; idx < _values.size(); idx++) {
      _values.at(idx) += other._values.at(idx);
    }

    return *this;
  }

  MatrixND operator+(const T &o) const {
    MatrixND ret(*this);

    for (uint64_t idx = 0; idx < ret._values.size(); idx++) {
      ret._values.at(idx) += o;
    }

    return ret;
  }

  MatrixND &operator+=(const T &o) {
    for (uint64_t idx = 0; idx < _values.size(); idx++) {
      _values.at(idx) += o;
    }

    return *this;
  }

  MatrixND operator-(const MatrixND &other) const {
    if (_sumSubDimensionCheck(other)) {
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

  MatrixND &operator-=(const MatrixND &other) {
    if (_sumSubDimensionCheck(other)) {
      std::ostringstream err;
      err << "Cannot substract matrices of dimensions: " << _printDimensions(_dimensions) << " and " << _printDimensions(other._dimensions);
      throw std::runtime_error(err.str());
    }

    for (uint64_t idx = 0; idx < _values.size(); idx++) {
      _values.at(idx) += other._values.at(idx);
    }

    return *this;
  }

  MatrixND operator-(const T &o) const {
    MatrixND ret(*this);

    for (uint64_t idx = 0; idx < ret._values.size(); idx++) {
      ret._values.at(idx) -= o;
    }

    return ret;
  }

  MatrixND &operator-=(const T &o) {
    for (uint64_t idx = 0; idx < _values.size(); idx++) {
      _values.at(idx) -= o;
    }

    return *this;
  }

  MatrixND operator*(const T &o) const {
    MatrixND ret(*this);

    for (uint64_t idx = 0; idx < ret._values.size(); idx++) {
      ret._values.at(idx) *= o;
    }

    return ret;
  }

  MatrixND &operator*=(const T &o) {
    for (uint64_t idx = 0; idx < _values.size(); idx++) {
      _values.at(idx) *= o;
    }

    return *this;
  }

  MatrixND operator/(const T &o) const {
    MatrixND ret(*this);

    for (uint64_t idx = 0; idx < ret._values.size(); idx++) {
      ret._values.at(idx) /= o;
    }

    return ret;
  }

  MatrixND &operator/=(const T &o) {
    for (uint64_t idx = 0; idx < _values.size(); idx++) {
      _values.at(idx) /= o;
    }

    return *this;
  }

  const std::vector<uint64_t> &getDimensions() const {
    return _dimensions;
  }

protected:
  const std::vector<uint64_t> _dimensions;

  std::vector<T> _values;

  bool _sumSubDimensionCheck(const MatrixND &other) const {
    return _dimensions != other._dimensions;
  }

  uint64_t _internalSize() const {
    uint64_t ret = 1;

    for (uint64_t idx = 0; idx < _dimensions.size(); idx++) {
      ret *= _dimensions[idx];
    }

    return ret;
  }

  uint64_t _indexOf(const std::vector<uint64_t> &coords) const {
    if (coords.size() != _dimensions.size()) {
      // NOTE: This shall never happen as we assure recourse creation is correct initialization
      // and this method is used in a protected context
      throw std::runtime_error("Invalid coordinates size");
    }

    for (uint64_t idx = 0; idx < coords.size(); idx++) {
      if (coords.at(idx) >= _dimensions.at(idx)) {
        std::ostringstream err;
        err << "Index {";
        for (const uint64_t c : coords) err << c << ", ";
        err << "} is out of range";
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

template<typename T = double>
class Matrix2D : public MatrixND<T> {
public:
  using MatrixND<T>::_indexOf;
  using MatrixND<T>::_printDimensions;

  Matrix2D() = delete;

  Matrix2D(uint64_t rows, uint64_t columns) : MatrixND<T>({rows, columns}), __rows(rows), __columns(columns) {}

  Matrix2D(uint64_t rows, uint64_t columns, const T &initial_value) : MatrixND<T>({rows, columns}, initial_value), __rows(rows), __columns(columns) {}

  Matrix2D(const Matrix2D &other) : MatrixND<T>(other), __rows(other.__rows), __columns(other.__columns) {}

  Matrix2D(Matrix2D &&other) : MatrixND<T>(std::move(other)), __rows(other.__rows), __columns(other.__columns) {}

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

  ~Matrix2D() = default;

  T &at(uint64_t row, uint64_t column) {
    return this->_values.at(_indexOf({row, column}));
  }

  const T &at(uint64_t row, uint64_t column) const {
    return this->_values.at(_indexOf({row, column}));
  }

  Matrix2D operator+(const Matrix2D &other) const {
    auto ret = MatrixND<T>::operator+(other);
    return Matrix2D::fromMatrixND(ret);
  }

  Matrix2D operator+(const T &o) const {
    Matrix2D ret(*this);
    ret += o;
    return ret;
  }

  Matrix2D operator-(const Matrix2D &other) const {
    auto ret = MatrixND<T>::operator-(other);
    return Matrix2D::fromMatrixND(ret);
  }

  Matrix2D operator-(const T &o) const {
    Matrix2D ret(*this);
    ret -= o;
    return ret;
  }

  Matrix2D operator*(const T &o) const {
    Matrix2D ret(*this);
    ret *= o;
    return ret;
  }

  Matrix2D operator/(const T &o) const {
    Matrix2D ret(*this);
    ret /= o;
    return ret;
  }

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

  uint64_t getRowAmount() const {
    return __rows;
  }

  uint64_t getColumnAmount() const {
    return __columns;
  }

  std::size_t getSize() const {
    return this->_values.size();
  }

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
  // Matrix2D(const MatrixND<T> &other) : MatrixND<T>(other), __rows(this->_dimensions.at(0)), __columns(this->_dimensions.at(1)) {}

  // Matrix2D(MatrixND<T> &&other) : MatrixND<T>(std::move(other)), __rows(this->_dimensions.at(0)), __columns(this->_dimensions.at(1)) {}

  const uint64_t __rows;
  const uint64_t __columns;
};

#endif /* __MATRIX_HPP__ */