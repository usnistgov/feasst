#ifndef FEASST_MATH_TABLE_H_
#define FEASST_MATH_TABLE_H_

#include <vector>
#include "utils/include/arguments.h"

namespace feasst {

/**
  Table values are assumed to be equally spaced and vary from 0 to 1.
 */
class Table {
 public:
  Table() {}

  /// Return the minimum of all elements.
  virtual double minimum() const = 0;

  /// Return the maximum of all elements.
  virtual double maximum() const = 0;

  /// Write to file.
  virtual void write(const std::string file_name) const;

  virtual ~Table() {}
};

/**
  This is a one-dimensional implementation of a table.
 */
class Table1D : public Table {
 public:
  /**
    args:
    - num: number of values (default: 1).
    - default_value: initialize data to this constant (default: 0).
   */
  Table1D(const argtype& args = argtype());

  /// Return the number of values.
  int num() const { return static_cast<int>(data_.size()); }

  /// For a given dimension, return the value of a bin.
  double bin_to_value(const int bin) const {
    return bin_spacing_*bin; }

  /// The inverse of above.
  int value_to_nearest_bin(const double value) const;

  /// Set data.
  void set_data(const int dim0, const double value) { data_[dim0] = value; }

  /// Return the data.
  const std::vector<double>& data() const { return data_; }

  /// Add the values of the given table.
  void add(const Table1D& table);

  /// Return linear interpolation of data given normalized values.
  double linear_interpolation(const double value0) const;

  /// Return the Newton-Gregory forward difference interpolation.
  /// See Allen and Tildesley, 5.2.2, Booth 1972
  double forward_difference_interpolation(const double value0) const;

  double minimum() const override;
  double maximum() const override;

  void serialize(std::ostream& ostr) const;
  explicit Table1D(std::istream& istr);

  // HWH python interface cannot handle stringstreams with serialization.
  std::string serialize() {
    std::stringstream ss;
    serialize(ss);
    return ss.str();
  }
  Table1D deserialize(const std::string str) {
    std::stringstream ss(str);
    return Table1D(ss);
  }

  virtual ~Table1D() {}

 private:
  std::vector<double> data_;
  double bin_spacing_;
  void calc_d_();
};

inline std::shared_ptr<Table1D> MakeTable1D(const argtype& args = argtype()) {
  return std::make_shared<Table1D>(args);
}

typedef std::vector<std::vector<double> > vec2;

/**
  This is a two-dimensional implementation of a table.
 */
class Table2D : public Table {
 public:
  /**
    args:
    - num0: number of values in first dimension (default: 1).
    - num1: number of values in second dimension (default: 1).
    - default_value: initialize data to this constant (default: 0).
   */
  Table2D(const argtype& args = argtype());

  /// Return the number of values in the first dimension
  int num0() const { return static_cast<int>(data_.size()); }

  /// Return the number of values in the second dimension
  int num1() const { return static_cast<int>(data_[0].size()); }

  /// Return the number of values in a given dimension.
  int num(const int dim) const;

  /// For a given dimension, return the value of a bin.
  double bin_to_value(const int dim, const int bin) const {
    return bin_spacing_[dim]*bin; }

  /// The inverse of above.
  int value_to_nearest_bin(const int dim, const double value) const;

  /// Set data.
  void set_data(const int dim0, const int dim1, const double value) {
    data_[dim0][dim1] = value; }

  /// Return the data.
  const vec2& data() const { return data_; }

  /// Add the values of the given table.
  void add(const Table2D& table);

  /// Return linear interpolation of data given normalized values for each
  /// dimension that range from 0 to 1, inclusive.
  double linear_interpolation(const double value0,
    const double value1) const;

  double minimum() const override;
  double maximum() const override;

  void serialize(std::ostream& ostr) const;
  explicit Table2D(std::istream& istr);

  // HWH python interface cannot handle stringstreams with serialization.
  std::string serialize() {
    std::stringstream ss;
    serialize(ss);
    return ss.str();
  }
  Table2D deserialize(const std::string str) {
    std::stringstream ss(str);
    return Table2D(ss);
  }

  virtual ~Table2D() {}

 private:
  vec2 data_;
  std::vector<double> bin_spacing_;
  void calc_d_();
};

inline std::shared_ptr<Table2D> MakeTable2D(const argtype& args = argtype()) {
  return std::make_shared<Table2D>(args);
}

typedef std::vector<std::vector<std::vector<double> > > vec3;

/**
  This is a three-dimensional implementation of a table.
 */
class Table3D : public Table {
 public:
  /**
    args:
    - num0: number of values in first dimension (default: 1).
    - num1: number of values in second dimension (default: 1).
    - num2: number of values in third dimension (default: 1).
    - default_value: initialize data to this constant (default: 0).
   */
  Table3D(const argtype& args = argtype());

  /// Return the number of values in the first dimension
  int num0() const { return static_cast<int>(data_.size()); }

  /// Return the number of values in the second dimension
  int num1() const { return static_cast<int>(data_[0].size()); }

  /// Return the number of values in the third dimension
  int num2() const { return static_cast<int>(data_[0][0].size()); }

  /// Return the number of values in a given dimension.
  int num(const int dim) const;

  /// For a given dimension, return the value of a bin.
  double bin_to_value(const int dim, const int bin) const {
    return bin_spacing_[dim]*bin; }

  /// The inverse of above.
  int value_to_nearest_bin(const int dim, const double value) const;

  /// Set data.
  void set_data(const int dim0, const int dim1, const int dim2,
    const double value) { data_[dim0][dim1][dim2] = value; }

  /// Return the data.
  const vec3& data() const { return data_; }

  /// Add the values of the given table.
  void add(const Table3D& table);

  /// Return linear interpolation of data given normalized values for each
  /// dimension that range from 0 to 1, inclusive.
  double linear_interpolation(const double value0,
    const double value1,
    const double value2) const;

  double minimum() const override;
  double maximum() const override;

  /// Write to file.
  void write(const std::string file_name) const override;

  /// Read from file.
  explicit Table3D(const std::string file_name);

  void serialize(std::ostream& ostr) const;
  explicit Table3D(std::istream& istr);

  // HWH python interface cannot handle stringstreams with serialization.
  std::string serialize() const {
    std::stringstream ss;
    serialize(ss);
    return ss.str();
  }
  Table3D deserialize(const std::string str) {
    std::stringstream ss(str);
    return Table3D(ss);
  }

  virtual ~Table3D() {}

 private:
  vec3 data_;
  std::vector<double> bin_spacing_;
  void calc_d_();
};

inline std::shared_ptr<Table3D> MakeTable3D(const argtype& args = argtype()) {
  return std::make_shared<Table3D>(args);
}

inline std::shared_ptr<Table3D> MakeTable3D(const std::string file_name) {
  return std::make_shared<Table3D>(file_name);
}

}  // namespace feasst

#endif  // FEASST_MATH_TABLE_H_
