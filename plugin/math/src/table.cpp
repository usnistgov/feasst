#include <string>
#include <fstream>
#include "utils/include/arguments.h"
#include "utils/include/utils.h"  // resize and fill
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "math/include/table.h"
#include "math/include/constants.h"

namespace feasst {

void Table::write(const std::string file_name) const {
  FATAL("not implemented");
}

double Table::calc_bin_spacing(const int num) {
  if (num > 1) {
    return 1./static_cast<double>(num - 1);
  } else {
    return 0;
  }
}

int Table::value_to_lowest_bin(const int dim, const double value) const {
  FATAL("Not Implemented.");
}

double Table::bin_to_value(const int dim, const int bin) const {
  FATAL("Not Implemented.");
}

void Table1D::calc_d_() {
  bin_spacing_ = calc_bin_spacing(num());
}

Table1D::Table1D(argtype args) : Table1D(&args) { feasst_check_all_used(args); }
Table1D::Table1D(argtype * args) : Table() {
  const int num = integer("num", args, 1);
  ASSERT(num > 0, "num:" << num << " must be > 0");
  data_.resize(num, dble("default_value", args, 0.));
  calc_d_();
}
Table1D::~Table1D() {}

double table_xd_(const double value0, const double d0, const int n0, int * i0, int * i02) {
  ASSERT(value0 >= 0 && value0 <= 1,
    "value: " << MAX_PRECISION << value0 << " out of bounds."
    << " If periodicity is disabled, place a hard wall at the boundaries.");
  TRACE("num0 " << n0);
  TRACE("value0 " << value0);
  // Numerical precision issue: If value is exactly 1, then i02 is exactly n0.
  if (value0 == 1.) {
    //*i02 = *i0;
    *i0 = n0 - 2;
    *i02 = n0 - 1;
    return 1.;
  }
  *i0 = value0*(n0 - 1);
  TRACE("i0 " << *i0);
  *i02 = *i0 + 1;
  TRACE("i02 " << *i02);
  const double v0 = *i0 * d0, vv0 = v0 + d0;
  TRACE("v0 " << v0);
  TRACE("vv0 " << vv0);
  TRACE("d0 " << d0);
  const double xd0 = (value0 - v0) / d0;
  TRACE("xd0 " << xd0);
  if (xd0 < 0) {
    if (-xd0 < 1e-14) {
      return 0;
    } else if (xd0 - 1. < 1e-14) {
      return 1.;
    }
  }
  ASSERT(xd0 >= 0 && xd0 <= 1, "xd0 " << xd0);
  return xd0;
}

double Table1D::c00_(const double xd0, const int i0, const int i02) const {
  const double c00 = data_[i0] * (1-xd0) + xd0*data_[i02];
  TRACE("c000 " << data_[i0]);
  TRACE("c100 " << data_[i02]);
  TRACE("c00 " << c00);
  return c00;
}

double Table1D::linear_interpolation(const double value0) const {
  int i0, i02;
  const double xd0 = table_xd_(value0, bin_spacing_, num(), &i0, &i02);
  return c00_(xd0, i0, i02);
}

double Table1D::forward_difference_interpolation(const double value0) const {
  const double sds = value0/bin_spacing_;
  TRACE("sds " << sds);
  const int k = static_cast<int>(sds);
  ASSERT(k < num(), "k: " << k << " beyond num: " << num());
  TRACE("k " << k);
  const double xi = sds - k;
  TRACE("xi " << xi);
  const double vk = data_[k];
  TRACE("vk " << vk);
  double vk1 = vk;
  if (k + 1 < num()) {
    vk1 = data_[k + 1];
  }
  TRACE("vk1 " << vk1);
  double vk2 = vk1;
  if (k + 2 < num()) {
    vk2 = data_[k + 2];
  }
  TRACE("vk2 " << vk2);
  const double t1 = vk + (vk1 - vk) * xi;
  TRACE("t1 " << t1);
  const double t2 = vk1 + (vk2 - vk1) * (xi - 1.);
  TRACE("t2 " << t2);
  const double rtrn_val = t1 + (t2 - t1)*xi*0.5;
  TRACE("rtrn_val: " << rtrn_val);
  return rtrn_val;
}

void Table1D::serialize(std::ostream& ostr) const {
  feasst_serialize_version(9075, ostr);
  feasst_serialize(data_, ostr);
}

Table1D::Table1D(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 9075, "version: " << version);
  feasst_deserialize(&data_, istr);
  calc_d_();
}

double Table1D::minimum() const { return feasst::minimum(data_); }
double Table1D::maximum() const { return feasst::maximum(data_); }

int Table1D::value_to_nearest_bin(const double value) const {
  return feasst::round(value*(num() - 1));
}

int Table1D::value_to_lowest_bin(const int dim, const double value) const {
  ASSERT(dim == 0, "dim:" << dim);
  int bin = static_cast<int>(value*(num() - 1));
  if (bin == num() - 1) {
    bin -= 1;
  } else if (bin == -1) {
    bin = 0;
  }
  TRACE("dim " << dim << " value " << value << " bin " << bin << " num " << num());
  return bin;
}

int Table1D::value_lowest_bin(const double value) const {
  return int(value*(num() - 1));
}

double Table1D::bin_to_value(const int dim, const int bin) const {
  ASSERT(dim == 0, "dim:" << dim);
  return bin_spacing_*bin;
}

void Table1D::add(const Table1D& table) { feasst::add(table.data_, &data_); }

void Table2D::calc_d_() {
  bin_spacing_ = std::vector<double>({
    calc_bin_spacing(num0()),
    calc_bin_spacing(num1())});
}

Table2D::Table2D(argtype args) : Table2D(&args) { feasst_check_all_used(args); }
Table2D::Table2D(argtype * args) : Table() {
  const int num0 = integer("num0", args, 1);
  const int num1 = integer("num1", args, 1);
  resize(num0, num1, &data_);
  fill(flt("default_value", args, 0.), &data_);
  calc_d_();
}

double Table2D::c00_(const double xd0, const double xd1, const int i0,
    const int i02, const int i1, const int i12) const {
  const double c00 = data_[i0][i1] *(1-xd0) + xd0*data_[i02][i1];
  const double c10 = data_[i0][i12]*(1-xd0) + xd0*data_[i02][i12];
  TRACE("c000 " << data_[i0][i1] << " c010 " << data_[i0][i12]);
  TRACE("c100 " << data_[i02][i1] << " c110 " << data_[i02][i12]);
  TRACE("c00 " << c00 << " c10 " << c10);
  const double c0 = c00 * (1-xd1) + xd1*c10;
  TRACE("c0 " << c0);
  return c0;
}

double Table2D::linear_interpolation(const double value0,
    const double value1) const {
  int i0, i02, i1, i12;
  const double xd0 = table_xd_(value0, bin_spacing_[0], num0(), &i0, &i02);
  const double xd1 = table_xd_(value1, bin_spacing_[1], num1(), &i1, &i12);
  return c00_(xd0, xd1, i0, i02, i1, i12);
}

void Table2D::serialize(std::ostream& ostr) const {
  feasst_serialize_version(2232, ostr);
  feasst_serialize(data_, ostr);
}

Table2D::Table2D(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2232, "version: " << version);
  feasst_deserialize(&data_, istr);
  calc_d_();
}

double Table2D::minimum() const { return feasst::minimum(data_); }
double Table2D::maximum() const { return feasst::maximum(data_); }

int Table2D::num(const int dim) const {
  if (dim == 0) {
    return num0();
  } else if (dim == 1) {
    return num1();
  } else {
    FATAL("dim: " << dim << " not recognized");
  }
}

int Table2D::value_to_nearest_bin(const int dim, const double value) const {
  return feasst::round(value*(num(dim) - 1));
}

int Table2D::value_to_lowest_bin(const int dim, const double value) const {
  int bin = static_cast<int>(value*(num(dim) - 1));
  if (bin == num(dim) - 1) {
    bin -= 1;
  } else if (bin == -1) {
    bin = 0;
  }
  TRACE("dim " << dim << " value " << value << " bin " << bin << " num " << num(dim));
  return bin;
}


void Table2D::add(const Table2D& table) { feasst::add(table.data_, &data_); }

void Table3D::calc_d_() {
  bin_spacing_ = std::vector<double>({
    calc_bin_spacing(num0()),
    calc_bin_spacing(num1()),
    calc_bin_spacing(num2())});
}

Table3D::Table3D(argtype args) : Table3D(&args) { feasst_check_all_used(args); }
Table3D::Table3D(argtype * args) : Table() {
  const int num0 = integer("num0", args, 1);
  const int num1 = integer("num1", args, 1);
  const int num2 = integer("num2", args, 1);
  resize(num0, num1, num2, &data_);
  fill(flt("default_value", args, 0.), &data_);
  calc_d_();
}

double Table3D::c00_(const double xd0, const double xd1, const double xd2,
    const int i0, const int i02, const int i1, const int i12, const int i2,
    const int i22) const {
  const double c00 = data_[i0][i1][i2] * (1-xd0) + xd0*data_[i02][i1][i2];
  const double c10 = data_[i0][i12][i2] *(1-xd0) + xd0*data_[i02][i12][i2];
  const double c01 = data_[i0][i1][i22] *(1-xd0) + xd0*data_[i02][i1][i22];
  const double c11 = data_[i0][i12][i22]*(1-xd0) + xd0*data_[i02][i12][i22];
  TRACE("c000 " << data_[i0][i1][i2] << " c010 " << data_[i0][i12][i2]
    << " c001 " << data_[i0][i1][i22] << " c011 " << data_[i0][i12][i22]);
  TRACE("c100 " << data_[i02][i1][i2] << " c110 " << data_[i02][i12][i2]
    << " c101 " << data_[i02][i1][i22] << " c111 " << data_[i02][i12][i22]);
  TRACE("c00 " << c00 << " c10 " << c10 << " c01 " << c01 << " c11 "
    << c11);

  const double c0 = c00 * (1-xd1) + xd1*c10;
  const double c1 = c01 * (1-xd1) + xd1*c11;
  TRACE("c0 " << c0 << " c1 " << c1);

  return c0*(1-xd2)+xd2*c1;
}

double Table3D::linear_interpolation(const double value0,
    const double value1, const double value2) const {
  int i0, i02, i1, i12, i2, i22;
  const double xd0 = table_xd_(value0, bin_spacing_[0], num0(), &i0, &i02);
  const double xd1 = table_xd_(value1, bin_spacing_[1], num1(), &i1, &i12);
  const double xd2 = table_xd_(value2, bin_spacing_[2], num2(), &i2, &i22);
  return c00_(xd0, xd1, xd2, i0, i02, i1, i12, i2, i22);
}

void Table3D::serialize(std::ostream& ostr) const {
  feasst_serialize_version(9826, ostr);
  feasst_serialize(data_, ostr);
}

Table3D::Table3D(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 9826 || version == 6867, "version: " << version);
  feasst_deserialize(&data_, istr);
  calc_d_();
}

double Table3D::minimum() const { return feasst::minimum(data_); }
double Table3D::maximum() const { return feasst::maximum(data_); }

int Table3D::num(const int dim) const {
  if (dim == 0) {
    return num0();
  } else if (dim == 1) {
    return num1();
  } else if (dim == 2) {
    return num2();
  } else {
    FATAL("dim: " << dim << " not recognized");
  }
}

int Table3D::value_to_nearest_bin(const int dim, const double value) const {
  return feasst::round(value*(num(dim) - 1));
}

int Table3D::value_to_lowest_bin(const int dim, const double value) const {
  int bin = static_cast<int>(value*(num(dim) - 1));
  if (bin == num(dim) - 1) {
    bin -= 1;
  } else if (bin == -1) {
    bin = 0;
  }
  return bin;
}

void Table3D::add(const Table3D& table) { feasst::add(table.data_, &data_); }

void Table3D::write(const std::string file_name) const {
  std::ofstream file(file_name);
  file << serialize();
}

Table3D::Table3D(const std::string file_name) {
  std::ifstream file(file_name);
  ASSERT(file.good(), "cannot find file " << file_name);
  std::string line;
  std::getline(file, line);
  *this = deserialize(line);
}

void Table4D::calc_d_() {
  bin_spacing_ = std::vector<double>({
    calc_bin_spacing(num0()),
    calc_bin_spacing(num1()),
    calc_bin_spacing(num2()),
    calc_bin_spacing(num3())});
}

Table4D::Table4D(argtype args) : Table4D(&args) { feasst_check_all_used(args); }
Table4D::Table4D(argtype * args) : Table() {
  const int num0 = integer("num0", args, 1);
  const int num1 = integer("num1", args, 1);
  const int num2 = integer("num2", args, 1);
  const int num3 = integer("num3", args, 1);
  resize(num0, num1, num2, num3, &data_);
  fill(flt("default_value", args, 0.), &data_);
  calc_d_();
}

double Table4D::linear_interpolation(const double value0,
    const double value1, const double value2, const double value3) const {
  int i0, i02, i1, i12, i2, i22, i3, i32;
  const double xd0 = table_xd_(value0, bin_spacing_[0], num0(), &i0, &i02);
  const double xd1 = table_xd_(value1, bin_spacing_[1], num1(), &i1, &i12);
  const double xd2 = table_xd_(value2, bin_spacing_[2], num2(), &i2, &i22);
  const double xd3 = table_xd_(value3, bin_spacing_[3], num3(), &i3, &i32);
  TRACE("size0 " << data_.size());
  TRACE("size1 " << data_[0].size());
  TRACE("size2 " << data_[0][0].size());
  TRACE("size3 " << data_[0][0][0].size());
  const double c000 = data_[i0][i1][i2][i3] * (1-xd0) +
                  xd0*data_[i02][i1][i2][i3];
  const double c100 = data_[i0][i12][i2][i3] *(1-xd0) +
                  xd0*data_[i02][i12][i2][i3];
  const double c010 = data_[i0][i1][i22][i3] *(1-xd0) +
                  xd0*data_[i02][i1][i22][i3];
  const double c110 = data_[i0][i12][i22][i3]*(1-xd0) +
                  xd0*data_[i02][i12][i22][i3];
  const double c001 = data_[i0][i1][i2][i32] * (1-xd0) +
                  xd0*data_[i02][i1][i2][i32];
  const double c101 = data_[i0][i12][i2][i32] *(1-xd0) +
                  xd0*data_[i02][i12][i2][i32];
  const double c011 = data_[i0][i1][i22][i32] *(1-xd0) +
                  xd0*data_[i02][i1][i22][i32];
  const double c111 = data_[i0][i12][i22][i32]*(1-xd0) +
                  xd0*data_[i02][i12][i22][i32];
  TRACE("c000 " << c000 << " c010 " << c010
    << " c001 " << c001 << " c011 " << c011
    << " c100 " << c100 << " c110 " << c110
    << " c101 " << c101 << " c111 " << c111);
  const double c00 = c000 * (1-xd1) + xd1*c100;
  const double c10 = c010 * (1-xd1) + xd1*c110;
  const double c01 = c001 * (1-xd1) + xd1*c101;
  const double c11 = c011 * (1-xd1) + xd1*c111;
  TRACE("c00 " << c00 << " c10 " << c10 << " c01 " << c01 << " c11 " << c11);
  const double c0 = c00*(1-xd2)+xd2*c10;
  const double c1 = c01*(1-xd2)+xd2*c11;
  TRACE("c0 " << c0 << " c1 " << c1);
  return c0*(1-xd3)+xd3*c1;
}

void Table4D::serialize(std::ostream& ostr) const {
  feasst_serialize_version(2707, ostr);
  feasst_serialize(data_, ostr);
}

Table4D::Table4D(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2707 || version == 6867, "version: " << version);
  feasst_deserialize(&data_, istr);
  calc_d_();
}

double Table4D::minimum() const { return feasst::minimum(data_); }
double Table4D::maximum() const { return feasst::maximum(data_); }

int Table4D::num(const int dim) const {
  if (dim == 0) {
    return num0();
  } else if (dim == 1) {
    return num1();
  } else if (dim == 2) {
    return num2();
  } else if (dim == 3) {
    return num3();
  } else {
    FATAL("dim: " << dim << " not recognized");
  }
}

int Table4D::value_to_nearest_bin(const int dim, const double value) const {
  return feasst::round(value*(num(dim) - 1));
}

void Table4D::add(const Table4D& table) { feasst::add(table.data_, &data_); }

void Table4D::write(const std::string file_name) const {
  std::ofstream file(file_name);
  file << serialize();
}

Table4D::Table4D(const std::string file_name) {
  std::ifstream file(file_name);
  ASSERT(file.good(), "cannot find file " << file_name);
  std::string line;
  std::getline(file, line);
  *this = deserialize(line);
}

void Table5D::calc_d_() {
  bin_spacing_ = std::vector<double>({
    calc_bin_spacing(num0()),
    calc_bin_spacing(num1()),
    calc_bin_spacing(num2()),
    calc_bin_spacing(num3()),
    calc_bin_spacing(num4())});
}

Table5D::Table5D(argtype args) : Table5D(&args) { feasst_check_all_used(args); }
Table5D::Table5D(argtype * args) : Table() {
  const int num0 = integer("num0", args, 1);
  const int num1 = integer("num1", args, 1);
  const int num2 = integer("num2", args, 1);
  const int num3 = integer("num3", args, 1);
  const int num4 = integer("num4", args, 1);
  resize(num0, num1, num2, num3, num4, &data_);
  fill(flt("default_value", args, 0.), &data_);
  calc_d_();
}

double Table5D::c00_(const double xd0, const double xd1, const double xd2,
    const double xd3, const double xd4, const int i0, const int i02,
    const int i1, const int i12, const int i2, const int i22, const int i3,
    const int i32, const int i4, const int i42) const {
  const double c0000 = data_[i0][i1][i2][i3][i4] * (1-xd0) +
                   xd0*data_[i02][i1][i2][i3][i4];
  const double c1000 = data_[i0][i12][i2][i3][i4] *(1-xd0) +
                   xd0*data_[i02][i12][i2][i3][i4];
  const double c0100 = data_[i0][i1][i22][i3][i4] *(1-xd0) +
                   xd0*data_[i02][i1][i22][i3][i4];
  const double c1100 = data_[i0][i12][i22][i3][i4]*(1-xd0) +
                   xd0*data_[i02][i12][i22][i3][i4];
  const double c0010 = data_[i0][i1][i2][i32][i4] * (1-xd0) +
                   xd0*data_[i02][i1][i2][i32][i4];
  const double c1010 = data_[i0][i12][i2][i32][i4] *(1-xd0) +
                   xd0*data_[i02][i12][i2][i32][i4];
  const double c0110 = data_[i0][i1][i22][i32][i4] *(1-xd0) +
                   xd0*data_[i02][i1][i22][i32][i4];
  const double c1110 = data_[i0][i12][i22][i32][i4]*(1-xd0) +
                   xd0*data_[i02][i12][i22][i32][i4];
  const double c0001 = data_[i0][i1][i2][i3][i42] * (1-xd0) +
                   xd0*data_[i02][i1][i2][i3][i42];
  const double c1001 = data_[i0][i12][i2][i3][i42] *(1-xd0) +
                   xd0*data_[i02][i12][i2][i3][i42];
  const double c0101 = data_[i0][i1][i22][i3][i42] *(1-xd0) +
                   xd0*data_[i02][i1][i22][i3][i42];
  const double c1101 = data_[i0][i12][i22][i3][i42]*(1-xd0) +
                   xd0*data_[i02][i12][i22][i3][i42];
  const double c0011 = data_[i0][i1][i2][i32][i42] * (1-xd0) +
                   xd0*data_[i02][i1][i2][i32][i42];
  const double c1011 = data_[i0][i12][i2][i32][i42] *(1-xd0) +
                   xd0*data_[i02][i12][i2][i32][i42];
  const double c0111 = data_[i0][i1][i22][i32][i42] *(1-xd0) +
                   xd0*data_[i02][i1][i22][i32][i42];
  const double c1111 = data_[i0][i12][i22][i32][i42]*(1-xd0) +
                   xd0*data_[i02][i12][i22][i32][i42];
  TRACE("c0000 " << c0000 << " c0100 " << c0100
    << " c0010 " << c0010 << " c0110 " << c0110
    << " c1000 " << c1000 << " c1100 " << c1100
    << " c1010 " << c1010 << " c1110 " << c1110
    << " c0001 " << c0001 << " c0101 " << c0101
    << " c0011 " << c0011 << " c0111 " << c0111
    << " c1001 " << c1001 << " c1101 " << c1101
    << " c1011 " << c1011 << " c1111 " << c1111);
  const double c000 = c0000*(1-xd1) + xd1*c1000;
  const double c100 = c0100*(1-xd1) + xd1*c1100;
  const double c010 = c0010*(1-xd1) + xd1*c1010;
  const double c110 = c0110*(1-xd1) + xd1*c1110;
  const double c001 = c0001*(1-xd1) + xd1*c1001;
  const double c101 = c0101*(1-xd1) + xd1*c1101;
  const double c011 = c0011*(1-xd1) + xd1*c1011;
  const double c111 = c0111*(1-xd1) + xd1*c1111;
  TRACE("c000 " << c000 << " c100 " << c100 <<
       " c010 " << c010 << " c110 " << c110 <<
        "c001 " << c001 << " c101 " << c101 <<
       " c011 " << c011 << " c111 " << c111);
  const double c00 = c000*(1-xd2)+xd2*c100;
  const double c10 = c010*(1-xd2)+xd2*c110;
  const double c01 = c001*(1-xd2)+xd2*c101;
  const double c11 = c011*(1-xd2)+xd2*c111;
  TRACE("c00 " << c00 << " c10 " << c10 <<
        "c01 " << c01 << " c11 " << c11);
  const double c0 = c00*(1-xd3)+xd3*c10;
  const double c1 = c01*(1-xd3)+xd3*c11;
  return c0*(1-xd4)+xd4*c1;
}

double Table5D::linear_interpolation(const double value0, const double value1,
    const double value2, const double value3, const double value4) const {
  int i0, i02, i1, i12, i2, i22, i3, i32, i4, i42;
  const double xd0 = table_xd_(value0, bin_spacing_[0], num0(), &i0, &i02);
  const double xd1 = table_xd_(value1, bin_spacing_[1], num1(), &i1, &i12);
  const double xd2 = table_xd_(value2, bin_spacing_[2], num2(), &i2, &i22);
  const double xd3 = table_xd_(value3, bin_spacing_[3], num3(), &i3, &i32);
  const double xd4 = table_xd_(value4, bin_spacing_[4], num4(), &i4, &i42);
  TRACE("size0 " << data_.size());
  TRACE("size1 " << data_[0].size());
  TRACE("size2 " << data_[0][0].size());
  TRACE("size3 " << data_[0][0][0].size());
  TRACE("size4 " << data_[0][0][0][0].size());
  return c00_(xd0, xd1, xd2, xd3, xd4, i0, i02, i1, i12, i2, i22, i3, i32, i4, i42);
}

void Table5D::serialize(std::ostream& ostr) const {
  feasst_serialize_version(6268, ostr);
  feasst_serialize(data_, ostr);
}

Table5D::Table5D(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6268 || version == 6867, "version: " << version);
  feasst_deserialize(&data_, istr);
  calc_d_();
}

double Table5D::minimum() const { return feasst::minimum(data_); }
double Table5D::maximum() const { return feasst::maximum(data_); }

int Table5D::num(const int dim) const {
  if (dim == 0) {
    return num0();
  } else if (dim == 1) {
    return num1();
  } else if (dim == 2) {
    return num2();
  } else if (dim == 3) {
    return num3();
  } else if (dim == 4) {
    return num4();
  } else {
    FATAL("dim: " << dim << " not recognized");
  }
}

int Table5D::value_to_nearest_bin(const int dim, const double value) const {
  return feasst::round(value*(num(dim) - 1));
}

int Table5D::value_to_lowest_bin(const int dim, const double value) const {
  int bin = static_cast<int>(value*(num(dim) - 1));
  if (bin == num(dim) - 1) {
    bin -= 1;
  } else if (bin == -1) {
    bin = 0;
  }
  return bin;
}

void Table5D::add(const Table5D& table) { feasst::add(table.data_, &data_); }

void Table5D::write(const std::string file_name) const {
  std::ofstream file(file_name);
  file << serialize();
}

Table5D::Table5D(const std::string file_name) {
  std::ifstream file(file_name);
  ASSERT(file.good(), "cannot find file " << file_name);
  std::string line;
  std::getline(file, line);
  *this = deserialize(line);
}

void Table6D::calc_d_() {
  bin_spacing_ = std::vector<double>({
    calc_bin_spacing(num0()),
    calc_bin_spacing(num1()),
    calc_bin_spacing(num2()),
    calc_bin_spacing(num3()),
    calc_bin_spacing(num4()),
    calc_bin_spacing(num5())});
}

Table6D::Table6D(argtype args) : Table6D(&args) { feasst_check_all_used(args); }
Table6D::Table6D(argtype * args) : Table() {
  const int num0 = integer("num0", args, 1);
  const int num1 = integer("num1", args, 1);
  const int num2 = integer("num2", args, 1);
  const int num3 = integer("num3", args, 1);
  const int num4 = integer("num4", args, 1);
  const int num5 = integer("num5", args, 1);
  resize(num0, num1, num2, num3, num4, num5, &data_);
  fill(flt("default_value", args, 0.), &data_);
  calc_d_();
}

double Table6D::c00_(const double xd0, const double xd1, const double xd2,
    const double xd3, const double xd4, const double xd5, const int i0, const int i02,
    const int i1, const int i12, const int i2, const int i22, const int i3,
    const int i32, const int i4, const int i42, const int i5, const int i52) const {
  const double c00000 = data_[i0][i1][i2][i3][i4][i5]     *(1-xd0) +
                    xd0*data_[i02][i1][i2][i3][i4][i5];
  const double c10000 = data_[i0][i12][i2][i3][i4][i5]    *(1-xd0) +
                    xd0*data_[i02][i12][i2][i3][i4][i5];
  const double c01000 = data_[i0][i1][i22][i3][i4][i5]    *(1-xd0) +
                    xd0*data_[i02][i1][i22][i3][i4][i5];
  const double c11000 = data_[i0][i12][i22][i3][i4][i5]   *(1-xd0) +
                    xd0*data_[i02][i12][i22][i3][i4][i5];
  const double c00100 = data_[i0][i1][i2][i32][i4][i5]    *(1-xd0) +
                    xd0*data_[i02][i1][i2][i32][i4][i5];
  const double c10100 = data_[i0][i12][i2][i32][i4][i5]   *(1-xd0) +
                    xd0*data_[i02][i12][i2][i32][i4][i5];
  const double c01100 = data_[i0][i1][i22][i32][i4][i5]   *(1-xd0) +
                    xd0*data_[i02][i1][i22][i32][i4][i5];
  const double c11100 = data_[i0][i12][i22][i32][i4][i5]  *(1-xd0) +
                    xd0*data_[i02][i12][i22][i32][i4][i5];
  const double c00010 = data_[i0][i1][i2][i3][i42][i5]    *(1-xd0) +
                    xd0*data_[i02][i1][i2][i3][i42][i5];
  const double c10010 = data_[i0][i12][i2][i3][i42][i5]   *(1-xd0) +
                    xd0*data_[i02][i12][i2][i3][i42][i5];
  const double c01010 = data_[i0][i1][i22][i3][i42][i5]   *(1-xd0) +
                    xd0*data_[i02][i1][i22][i3][i42][i5];
  const double c11010 = data_[i0][i12][i22][i3][i42][i5]  *(1-xd0) +
                    xd0*data_[i02][i12][i22][i3][i42][i5];
  const double c00110 = data_[i0][i1][i2][i32][i42][i5]   *(1-xd0) +
                    xd0*data_[i02][i1][i2][i32][i42][i5];
  const double c10110 = data_[i0][i12][i2][i32][i42][i5]  *(1-xd0) +
                    xd0*data_[i02][i12][i2][i32][i42][i5];
  const double c01110 = data_[i0][i1][i22][i32][i42][i5]  *(1-xd0) +
                    xd0*data_[i02][i1][i22][i32][i42][i5];
  const double c11110 = data_[i0][i12][i22][i32][i42][i5] *(1-xd0) +
                    xd0*data_[i02][i12][i22][i32][i42][i5];
  const double c00001 = data_[i0][i1][i2][i3][i4][i52]    *(1-xd0) +
                    xd0*data_[i02][i1][i2][i3][i4][i52];
  const double c10001 = data_[i0][i12][i2][i3][i4][i52]   *(1-xd0) +
                    xd0*data_[i02][i12][i2][i3][i4][i52];
  const double c01001 = data_[i0][i1][i22][i3][i4][i52]   *(1-xd0) +
                    xd0*data_[i02][i1][i22][i3][i4][i52];
  const double c11001 = data_[i0][i12][i22][i3][i4][i52]  *(1-xd0) +
                    xd0*data_[i02][i12][i22][i3][i4][i52];
  const double c00101 = data_[i0][i1][i2][i32][i4][i52]   *(1-xd0) +
                    xd0*data_[i02][i1][i2][i32][i4][i52];
  const double c10101 = data_[i0][i12][i2][i32][i4][i52]  *(1-xd0) +
                    xd0*data_[i02][i12][i2][i32][i4][i52];
  const double c01101 = data_[i0][i1][i22][i32][i4][i52]  *(1-xd0) +
                    xd0*data_[i02][i1][i22][i32][i4][i52];
  const double c11101 = data_[i0][i12][i22][i32][i4][i52] *(1-xd0) +
                    xd0*data_[i02][i12][i22][i32][i4][i52];
  const double c00011 = data_[i0][i1][i2][i3][i42][i52]   *(1-xd0) +
                    xd0*data_[i02][i1][i2][i3][i42][i52];
  const double c10011 = data_[i0][i12][i2][i3][i42][i52]  *(1-xd0) +
                    xd0*data_[i02][i12][i2][i3][i42][i52];
  const double c01011 = data_[i0][i1][i22][i3][i42][i52]  *(1-xd0) +
                    xd0*data_[i02][i1][i22][i3][i42][i52];
  const double c11011 = data_[i0][i12][i22][i3][i42][i52] *(1-xd0) +
                    xd0*data_[i02][i12][i22][i3][i42][i52];
  const double c00111 = data_[i0][i1][i2][i32][i42][i52]  *(1-xd0) +
                    xd0*data_[i02][i1][i2][i32][i42][i52];
  const double c10111 = data_[i0][i12][i2][i32][i42][i52] *(1-xd0) +
                    xd0*data_[i02][i12][i2][i32][i42][i52];
  const double c01111 = data_[i0][i1][i22][i32][i42][i52] *(1-xd0) +
                    xd0*data_[i02][i1][i22][i32][i42][i52];
  const double c11111 = data_[i0][i12][i22][i32][i42][i52]*(1-xd0) +
                    xd0*data_[i02][i12][i22][i32][i42][i52];
  const double c0000 = c00000*(1-xd1) + xd1*c10000;
  const double c1000 = c01000*(1-xd1) + xd1*c11000;
  const double c0100 = c00100*(1-xd1) + xd1*c10100;
  const double c1100 = c01100*(1-xd1) + xd1*c11100;
  const double c0010 = c00010*(1-xd1) + xd1*c10010;
  const double c1010 = c01010*(1-xd1) + xd1*c11010;
  const double c0110 = c00110*(1-xd1) + xd1*c10110;
  const double c1110 = c01110*(1-xd1) + xd1*c11110;
  const double c0001 = c00001*(1-xd1) + xd1*c10001;
  const double c1001 = c01001*(1-xd1) + xd1*c11001;
  const double c0101 = c00101*(1-xd1) + xd1*c10101;
  const double c1101 = c01101*(1-xd1) + xd1*c11101;
  const double c0011 = c00011*(1-xd1) + xd1*c10011;
  const double c1011 = c01011*(1-xd1) + xd1*c11011;
  const double c0111 = c00111*(1-xd1) + xd1*c10111;
  const double c1111 = c01111*(1-xd1) + xd1*c11111;
  const double c000 = c0000*(1-xd2)+xd2*c1000;
  const double c100 = c0100*(1-xd2)+xd2*c1100;
  const double c010 = c0010*(1-xd2)+xd2*c1010;
  const double c110 = c0110*(1-xd2)+xd2*c1110;
  const double c001 = c0001*(1-xd2)+xd2*c1001;
  const double c101 = c0101*(1-xd2)+xd2*c1101;
  const double c011 = c0011*(1-xd2)+xd2*c1011;
  const double c111 = c0111*(1-xd2)+xd2*c1111;
  const double c00 = c000*(1-xd3)+xd3*c100;
  const double c10 = c010*(1-xd3)+xd3*c110;
  const double c01 = c001*(1-xd3)+xd3*c101;
  const double c11 = c011*(1-xd3)+xd3*c111;
  const double c0 = c00*(1-xd4)+xd4*c10;
  const double c1 = c01*(1-xd4)+xd4*c11;
  const double rtn = c0*(1-xd5)+xd5*c1;
  return rtn;
}

double Table6D::linear_interpolation(const double value0, const double value1,
    const double value2, const double value3, const double value4,
    const double value5) const {
  int i0, i02, i1, i12, i2, i22, i3, i32, i4, i42, i5, i52;
  const double xd0 = table_xd_(value0, bin_spacing_[0], num0(), &i0, &i02);
  const double xd1 = table_xd_(value1, bin_spacing_[1], num1(), &i1, &i12);
  const double xd2 = table_xd_(value2, bin_spacing_[2], num2(), &i2, &i22);
  const double xd3 = table_xd_(value3, bin_spacing_[3], num3(), &i3, &i32);
  const double xd4 = table_xd_(value4, bin_spacing_[4], num4(), &i4, &i42);
  const double xd5 = table_xd_(value5, bin_spacing_[5], num5(), &i5, &i52);
  TRACE("xd0 " << xd0 << " xd1 " << xd1 << " xd2 " << xd2 << " xd3 " << xd3 <<
       " xd4 " << xd4 << " xd5 " << xd5);
  TRACE("size0 " << data_.size());
  TRACE("size1 " << data_[0].size());
  TRACE("size2 " << data_[0][0].size());
  TRACE("size3 " << data_[0][0][0].size());
  TRACE("size4 " << data_[0][0][0][0].size());
  TRACE("size5 " << data_[0][0][0][0][0].size());
  return c00_(xd0, xd1, xd2, xd3, xd4, xd5, i0, i02, i1, i12, i2, i22, i3, i32, i4, i42, i5, i52);
//  if (std::isnan(rtn)) {
//    INFO("value0 " << value0 << "value1 " << value1 << "value2 " << value2 <<
//      " value3 " << value3 << " value4 " << value4 << " value5 " << value5);
//    INFO("i0 " << i0 << " i1 " << i1 << " i2 " << i2 << " i3 " << i3 << " i4 "
//      << i4 << " i5 " << i5);
//    INFO("i02 " << i02 << " i12 " << i12 << " i22 " << i22 << " i32 " << i32 <<
//      " i42 " << i42 << " i52 " << i52);
//    INFO("xd0 " << xd0 << " xd1 " << xd1 << " xd2 " << xd2 << " xd3 " << xd3 <<
//      " xd4 " << xd4 << " xd5 " << xd5);
//    INFO("size0 " << data_.size());
//    INFO("size1 " << data_[0].size());
//    INFO("size2 " << data_[0][0].size());
//    INFO("size3 " << data_[0][0][0].size());
//    INFO("size4 " << data_[0][0][0][0].size());
//    INFO("size5 " << data_[0][0][0][0][0].size());
//    INFO("c0 " << c0 << " c1 " << c1);
//    INFO("c01 " << c01 << " c11 " << c11);
//    INFO("c011 " << c011 << " c111 " << c111);
//    INFO("c0111 " << c0111 << " c1111 " << c1111);
//    INFO("c01111 " << c01111 << " c11111 " << c11111);
//    INFO(data_[i0][i12][i22][i32][i42][i52]);
//    INFO(data_[i02][i12][i22][i32][i42][i52]);
//    INFO(i02 << " " << i12 << " " << i22 << " " << i32 << " " << i42 << " "
//      << i52);
//    FATAL("fatal");
//  }
//  return rtn;
}

void Table6D::serialize(std::ostream& ostr) const {
  feasst_serialize_version(6867, ostr);
  feasst_serialize(data_, ostr);
}

Table6D::Table6D(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6867, "version: " << version);
  feasst_deserialize(&data_, istr);
  calc_d_();
}

double Table6D::minimum() const { return feasst::minimum(data_); }
double Table6D::maximum() const { return feasst::maximum(data_); }

int Table6D::num(const int dim) const {
  if (dim == 0) {
    return num0();
  } else if (dim == 1) {
    return num1();
  } else if (dim == 2) {
    return num2();
  } else if (dim == 3) {
    return num3();
  } else if (dim == 4) {
    return num4();
  } else if (dim == 5) {
    return num5();
  } else {
    FATAL("dim: " << dim << " not recognized");
  }
}

int Table6D::value_to_nearest_bin(const int dim, const double value) const {
  return feasst::round(value*(num(dim) - 1));
}

int Table6D::value_to_lowest_bin(const int dim, const double value) const {
  int bin = static_cast<int>(value*(num(dim) - 1));
  if (bin == num(dim) - 1) {
    bin -= 1;
  } else if (bin == -1) {
    bin = 0;
  }
  return bin;
}

void Table6D::add(const Table6D& table) { feasst::add(table.data_, &data_); }

void Table6D::write(const std::string file_name) const {
  std::ofstream file(file_name);
  file << serialize();
}

Table6D::Table6D(const std::string file_name) {
  std::ifstream file(file_name);
  ASSERT(file.good(), "cannot find file " << file_name);
  std::string line;
  std::getline(file, line);
  *this = deserialize(line);
}

std::string Table1D::serialize() const {
  std::stringstream ss;
  serialize(ss);
  return ss.str();
}
Table1D Table1D::deserialize(const std::string str) {
  std::stringstream ss(str);
  return Table1D(ss);
}
std::string Table2D::serialize() const {
  std::stringstream ss;
  serialize(ss);
  return ss.str();
}
Table2D Table2D::deserialize(const std::string str) {
  std::stringstream ss(str);
  return Table2D(ss);
}
std::string Table3D::serialize() const {
  std::stringstream ss;
  serialize(ss);
  return ss.str();
}
Table3D Table3D::deserialize(const std::string str) {
  std::stringstream ss(str);
  return Table3D(ss);
}
std::string Table4D::serialize() const {
  std::stringstream ss;
  serialize(ss);
  return ss.str();
}
Table4D Table4D::deserialize(const std::string str) {
  std::stringstream ss(str);
  return Table4D(ss);
}
std::string Table5D::serialize() const {
  std::stringstream ss;
  serialize(ss);
  return ss.str();
}
Table5D Table5D::deserialize(const std::string str) {
  std::stringstream ss(str);
  return Table5D(ss);
}
std::string Table6D::serialize() const {
  std::stringstream ss;
  serialize(ss);
  return ss.str();
}
Table6D Table6D::deserialize(const std::string str) {
  std::stringstream ss(str);
  return Table6D(ss);
}
}  // namespace feasst
