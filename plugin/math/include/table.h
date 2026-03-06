#ifndef FEASST_MATH_TABLE_H_
#define FEASST_MATH_TABLE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace feasst {

typedef std::map<std::string, std::string> argtype;

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

  /// Return the bin spacing for a given number of elements.
  double calc_bin_spacing(const int num);

  /// Write to file.
  virtual void write(const std::string file_name) const;

  virtual ~Table() {}
};

double table_xd_(const double value0, const double d0, const int n0, int * i0, int * i02);

/**
  One-dimensional implementation of a table.
 */
class Table1D : public Table {
 public:
  /**
    args:
    - num: number of values (default: 1).
    - default_value: initialize data to this constant (default: 0).
   */
  explicit Table1D(argtype args = argtype());
  explicit Table1D(argtype *args);

  /// Return the number of values.
  int num() const { return static_cast<int>(data_.size()); }

  /// Return the bin spacing.
  double bin_spacing() const { return bin_spacing_; }

  /// For a given dimension, return the value of a bin.
  double bin_to_value(const int bin) const {
    return bin_spacing_*bin; }

  /// The inverse of above.
  int value_to_nearest_bin(const double value) const;

  /// Return the bin just below the value.
  int value_to_lowest_bin(const double value) const;

  /// Set data.
  void set_data(const int dim0, const double value) { data_[dim0] = value; }

  /// Return the data.
  const std::vector<double>& data() const { return data_; }

  /// Return the data.
  double data(const int index) const { return data_[index]; }

  /// Add the values of the given table.
  void add(const Table1D& table);

  /// Return linear interpolation of data given normalized values.
  virtual double linear_interpolation(const double value0) const;

  /// Return the Newton-Gregory forward difference interpolation.
  /// See Allen and Tildesley, 5.2.2, Booth 1972
  virtual double forward_difference_interpolation(const double value0) const;

  double minimum() const override;
  double maximum() const override;

  void serialize(std::ostream& ostr) const;
  explicit Table1D(std::istream& istr);

  // HWH python interface cannot handle stringstreams with serialization.
  std::string serialize() const;
  Table1D deserialize(const std::string str);

  virtual ~Table1D();

 protected:
  double c00_(const double xd0, const int i0, const int i02) const;

 private:
  std::vector<double> data_;
  double bin_spacing_;
  void calc_d_();
};

inline std::shared_ptr<Table1D> MakeTable1D(argtype args = argtype()) {
  return std::make_shared<Table1D>(args); }

typedef std::vector<std::vector<float> > fvec2;

/**
  Two-dimensional implementation of a table.
 */
class Table2D : public Table {
 public:
  /**
    args:
    - num0: number of values in first dimension (default: 1).
    - num1: number of values in second dimension (default: 1).
    - default_value: initialize data to this constant (default: 0).
   */
  explicit Table2D(argtype args = argtype());
  explicit Table2D(argtype * args);

  /// Return the number of values in the first dimension
  int num0() const { return static_cast<int>(data_.size()); }

  /// Return the number of values in the second dimension
  int num1() const { return static_cast<int>(data_[0].size()); }

  /// Return the number of values in a given dimension.
  int num(const int dim) const;

  /// Return the bin spacing.
  double bin_spacing(const int dim) const { return bin_spacing_[dim]; }

  /// For a given dimension, return the value of a bin.
  double bin_to_value(const int dim, const int bin) const {
    return bin_spacing_[dim]*bin; }

  /// The inverse of above.
  int value_to_nearest_bin(const int dim, const double value) const;
  int value_to_lowest_bin(const int dim, const double value) const;

  /// Set data.
  void set_data(const int dim0, const int dim1, const double value) {
    data_[dim0][dim1] = value; }

  /// Return the data.
  const fvec2& data() const { return data_; }

  /// Add the values of the given table.
  void add(const Table2D& table);

  /// Return linear interpolation of data given normalized values for each
  /// dimension that range from 0 to 1, inclusive.
  virtual double linear_interpolation(const double value0,
    const double value1) const;

  double minimum() const override;
  double maximum() const override;

  void serialize(std::ostream& ostr) const;
  explicit Table2D(std::istream& istr);

  // HWH python interface cannot handle stringstreams with serialization.
  std::string serialize() const;
  Table2D deserialize(const std::string str);

  virtual ~Table2D() {}

 protected:
  double c00_(const double xd0, const double xd1, const int i0, const int i02,
    const int i1, const int i12) const;

 private:
  fvec2 data_;
  std::vector<double> bin_spacing_;
  void calc_d_();
};

inline std::shared_ptr<Table2D> MakeTable2D(argtype args = argtype()) {
  return std::make_shared<Table2D>(args); }

typedef std::vector<fvec2> fvec3;

/**
  Three-dimensional implementation of a table.
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
  explicit Table3D(argtype args = argtype());
  explicit Table3D(argtype * args);

  /// Return the number of values in the first dimension
  int num0() const { return static_cast<int>(data_.size()); }

  /// Return the number of values in the second dimension
  int num1() const { return static_cast<int>(data_[0].size()); }

  /// Return the number of values in the third dimension
  int num2() const { return static_cast<int>(data_[0][0].size()); }

  /// Return the number of values in a given dimension.
  int num(const int dim) const;

  /// Return the bin spacing.
  double bin_spacing(const int dim) const { return bin_spacing_[dim]; }

  /// For a given dimension, return the value of a bin.
  double bin_to_value(const int dim, const int bin) const {
    return bin_spacing_[dim]*bin; }

  /// The inverse of above.
  int value_to_nearest_bin(const int dim, const double value) const;

  /// Set data.
  void set_data(const int dim0, const int dim1, const int dim2,
    const double value) { data_[dim0][dim1][dim2] = value; }

  /// Return the data.
  const fvec3& data() const { return data_; }

  /// Add the values of the given table.
  void add(const Table3D& table);

  /// Return linear interpolation of data given normalized values for each
  /// dimension that range from 0 to 1, inclusive.
  virtual double linear_interpolation(const double value0,
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
  std::string serialize() const;
  Table3D deserialize(const std::string str);

  virtual ~Table3D() {}

 protected:
  double c00_(const double xd0, const double xd1, const double xd2,
    const int i0, const int i02, const int i1, const int i12, const int i2,
    const int i22) const;

 private:
  fvec3 data_;
  std::vector<double> bin_spacing_;
  void calc_d_();
};

inline std::shared_ptr<Table3D> MakeTable3D(argtype args = argtype()) {
  return std::make_shared<Table3D>(args); }

inline std::shared_ptr<Table3D> MakeTable3D(const std::string file_name) {
  return std::make_shared<Table3D>(file_name);
}

typedef std::vector<fvec3> fvec4;

/**
  Four-dimensional implementation of a table.
 */
class Table4D : public Table {
 public:
  /**
    args:
    - num0: number of values in first dimension (default: 1).
    - num1: number of values in second dimension (default: 1).
    - num2: number of values in third dimension (default: 1).
    - num3: number of values in fourth dimension (default: 1).
    - default_value: initialize data to this constant (default: 0).
   */
  explicit Table4D(argtype args = argtype());
  explicit Table4D(argtype * args);

  /// Return the number of values in the first dimension
  int num0() const { return static_cast<int>(data_.size()); }

  /// Return the number of values in the second dimension
  int num1() const { return static_cast<int>(data_[0].size()); }

  /// Return the number of values in the third dimension
  int num2() const { return static_cast<int>(data_[0][0].size()); }

  /// Return the number of values in the fourth dimension
  int num3() const { return static_cast<int>(data_[0][0][0].size()); }

  /// Return the number of values in a given dimension.
  int num(const int dim) const;

  /// Return the bin spacing.
  double bin_spacing(const int dim) const { return bin_spacing_[dim]; }

  /// For a given dimension, return the value of a bin.
  double bin_to_value(const int dim, const int bin) const {
    return bin_spacing_[dim]*bin; }

  /// The inverse of above.
  int value_to_nearest_bin(const int dim, const double value) const;

  /// Set data.
  void set_data(const int dim0, const int dim1, const int dim2, const int dim3,
    const double value) { data_[dim0][dim1][dim2][dim3] = value; }

  /// Return the data.
  const fvec4& data() const { return data_; }

  /// Add the values of the given table.
  void add(const Table4D& table);

  /// Return linear interpolation of data given normalized values for each
  /// dimension that range from 0 to 1, inclusive.
  double linear_interpolation(const double value0,
    const double value1,
    const double value2,
    const double value3) const;

  double minimum() const override;
  double maximum() const override;

  /// Write to file.
  void write(const std::string file_name) const override;

  /// Read from file.
  explicit Table4D(const std::string file_name);

  void serialize(std::ostream& ostr) const;
  explicit Table4D(std::istream& istr);

  // HWH python interface cannot handle stringstreams with serialization.
  std::string serialize() const;
  Table4D deserialize(const std::string str);

  virtual ~Table4D() {}

 private:
  fvec4 data_;
  std::vector<double> bin_spacing_;
  void calc_d_();
};

inline std::shared_ptr<Table4D> MakeTable4D(argtype args = argtype()) {
  return std::make_shared<Table4D>(args); }

inline std::shared_ptr<Table4D> MakeTable4D(const std::string file_name) {
  return std::make_shared<Table4D>(file_name);
}

typedef std::vector<fvec4> fvec5;

/**
  Five-dimensional implementation of a table.
 */
class Table5D : public Table {
 public:
  /**
    args:
    - num0: number of values in first dimension (default: 1).
    - num1: number of values in second dimension (default: 1).
    - num2: number of values in third dimension (default: 1).
    - num3: number of values in fourth dimension (default: 1).
    - num4: number of values in fourth dimension (default: 1).
    - default_value: initialize data to this constant (default: 0).
   */
  explicit Table5D(argtype args = argtype());
  explicit Table5D(argtype * args);

  /// Return the number of values in the first dimension
  int num0() const { return static_cast<int>(data_.size()); }

  /// Return the number of values in the second dimension
  int num1() const { return static_cast<int>(data_[0].size()); }

  /// Return the number of values in the third dimension
  int num2() const { return static_cast<int>(data_[0][0].size()); }

  /// Return the number of values in the fourth dimension
  int num3() const { return static_cast<int>(data_[0][0][0].size()); }

  /// Return the number of values in the fifth dimension
  int num4() const { return static_cast<int>(data_[0][0][0][0].size()); }

  /// Return the number of values in a given dimension.
  int num(const int dim) const;

  /// Return the bin spacing.
  double bin_spacing(const int dim) const { return bin_spacing_[dim]; }

  /// For a given dimension, return the value of a bin.
  double bin_to_value(const int dim, const int bin) const {
    return bin_spacing_[dim]*bin; }

  /// The inverse of above.
  int value_to_nearest_bin(const int dim, const double value) const;
  int value_to_lowest_bin(const int dim, const double value) const;

  /// Set data.
  void set_data(const int dim0, const int dim1, const int dim2, const int dim3,
      const int dim4, const double value) {
    data_[dim0][dim1][dim2][dim3][dim4] = value; }

  /// Return the data.
  const fvec5& data() const { return data_; }

  /// Add the values of the given table.
  void add(const Table5D& table);

  /// Return linear interpolation of data given normalized values for each
  /// dimension that range from 0 to 1, inclusive.
  virtual double linear_interpolation(const double value0,
    const double value1,
    const double value2,
    const double value3,
    const double value4) const;

  double minimum() const override;
  double maximum() const override;

  /// Write to file.
  void write(const std::string file_name) const override;

  /// Read from file.
  explicit Table5D(const std::string file_name);

  void serialize(std::ostream& ostr) const;
  explicit Table5D(std::istream& istr);

  // HWH python interface cannot handle stringstreams with serialization.
  std::string serialize() const;
  Table5D deserialize(const std::string str);

  virtual ~Table5D() {}
 protected:
  double c00_(const double xd0, const double xd1, const double xd2,
    const double xd3, const double xd4, const int i0, const int i02,
    const int i1, const int i12, const int i2, const int i22, const int i3,
    const int i32, const int i4, const int i42) const;

 private:
  fvec5 data_;
  std::vector<double> bin_spacing_;
  void calc_d_();
};

inline std::shared_ptr<Table5D> MakeTable5D(argtype args = argtype()) {
  return std::make_shared<Table5D>(args); }

inline std::shared_ptr<Table5D> MakeTable5D(const std::string file_name) {
  return std::make_shared<Table5D>(file_name);
}

typedef std::vector<fvec5> fvec6;

/**
  Six-dimensional implementation of a table.
 */
class Table6D : public Table {
 public:
  /**
    args:
    - num0: number of values in first dimension (default: 1).
    - num1: number of values in second dimension (default: 1).
    - num2: number of values in third dimension (default: 1).
    - num3: number of values in fourth dimension (default: 1).
    - num4: number of values in fourth dimension (default: 1).
    - num5: number of values in fourth dimension (default: 1).
    - default_value: initialize data to this constant (default: 0).
   */
  explicit Table6D(argtype args = argtype());
  explicit Table6D(argtype * args);

  /// Return the number of values in the first dimension
  int num0() const { return static_cast<int>(data_.size()); }

  /// Return the number of values in the second dimension
  int num1() const { return static_cast<int>(data_[0].size()); }

  /// Return the number of values in the third dimension
  int num2() const { return static_cast<int>(data_[0][0].size()); }

  /// Return the number of values in the fourth dimension
  int num3() const { return static_cast<int>(data_[0][0][0].size()); }

  /// Return the number of values in the fifth dimension
  int num4() const { return static_cast<int>(data_[0][0][0][0].size()); }

  /// Return the number of values in the sixth dimension
  int num5() const { return static_cast<int>(data_[0][0][0][0][0].size()); }

  /// Return the number of values in a given dimension.
  int num(const int dim) const;

  /// Return the bin spacing.
  double bin_spacing(const int dim) const { return bin_spacing_[dim]; }

  /// For a given dimension, return the value of a bin.
  double bin_to_value(const int dim, const int bin) const {
    return bin_spacing_[dim]*bin; }

  /// The inverse of above.
  int value_to_nearest_bin(const int dim, const double value) const;

  /// Set data.
  void set_data(const int dim0, const int dim1, const int dim2, const int dim3,
      const int dim4, const int dim5, const double value) {
    data_[dim0][dim1][dim2][dim3][dim4][dim5] = value; }

  /// Return the data.
  const fvec6& data() const { return data_; }

  /// Add the values of the given table.
  void add(const Table6D& table);

  /// Return linear interpolation of data given normalized values for each
  /// dimension that range from 0 to 1, inclusive.
  virtual double linear_interpolation(const double value0,
    const double value1,
    const double value2,
    const double value3,
    const double value4,
    const double value5) const;

  double minimum() const override;
  double maximum() const override;

  /// Write to file.
  void write(const std::string file_name) const override;

  /// Read from file.
  explicit Table6D(const std::string file_name);

  void serialize(std::ostream& ostr) const;
  explicit Table6D(std::istream& istr);

  // HWH python interface cannot handle stringstreams with serialization.
  std::string serialize() const;
  Table6D deserialize(const std::string str);

  virtual ~Table6D() {}

 protected:
  double c00_(const double xd0, const double xd1, const double xd2,
    const double xd3, const double xd4, const double xd5, const int i0, const int i02,
    const int i1, const int i12, const int i2, const int i22, const int i3,
    const int i32, const int i4, const int i42, const int i5, const int i52) const;

 private:
  fvec6 data_;
  std::vector<double> bin_spacing_;
  void calc_d_();
};

inline std::shared_ptr<Table6D> MakeTable6D(argtype args = argtype()) {
  return std::make_shared<Table6D>(args); }

inline std::shared_ptr<Table6D> MakeTable6D(const std::string file_name) {
  return std::make_shared<Table6D>(file_name);
}

}  // namespace feasst

#endif  // FEASST_MATH_TABLE_H_
