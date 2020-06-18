#ifndef FEASST_MATH_TABLE_H_
#define FEASST_MATH_TABLE_H_

#include <vector>
#include "utils/include/arguments.h"

namespace feasst {

typedef std::vector<std::vector<std::vector<double> > > vec3;

/**
  Table values are assumed to be equally spaced and vary from 0 to 1.
  This is a three-dimensional implementation of a table.
 */
class Table3D {
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

  /// For a given dimension, return the value of a bin.
  double bin_to_value(const int dim, const int bin) const {
    return bin_spacing_[dim]*bin; }

  /// Set data.
  void set_data(const int dim0, const int dim1, const int dim2,
    const double value) { data_[dim0][dim1][dim2] = value; }

  /// Return the data.
  const vec3& data() const { return data_; }

  /// Return linear interpolation of data given normalized values for each
  /// dimension that range from 0 to 1, inclusive.
  double linear_interpolation(const double value0,
    const double value1,
    const double value2) const;

  void serialize(std::ostream& ostr) const;
  explicit Table3D(std::istream& istr);

  // HWH python interface cannot handle stringstreams with serialization.
  std::string serialize() {
    std::stringstream ss;
    serialize(ss);
    return ss.str();
  }
  Table3D deserialize(const std::string str) {
    std::stringstream ss(str);
    return Table3D(ss);
  }

 private:
  vec3 data_;
  std::vector<double> bin_spacing_;

  void calc_d_();
};

inline std::shared_ptr<Table3D> MakeTable3D(const argtype& args = argtype()) {
  return std::make_shared<Table3D>(args);
}

}  // namespace feasst

#endif  // FEASST_MATH_TABLE_H_
