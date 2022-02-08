#include <string>
#include <fstream>
#include "utils/include/utils.h"  // resize and fill
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "math/include/table.h"
#include "math/include/constants.h"

namespace feasst {

void Table::write(const std::string file_name) const {
  FATAL("not implemented");
}

double Table::bin_spacing(const int num) {
  if (num > 1) {
    return 1./static_cast<double>(num - 1);
  } else {
    return 0;
  }
}

void Table1D::calc_d_() {
  bin_spacing_ = bin_spacing(num());
}

Table1D::Table1D(argtype args) : Table1D(&args) { check_all_used(args); }
Table1D::Table1D(argtype * args) : Table() {
  const int num = integer("num", args, 1);
  data_.resize(num, dble("default_value", args, 0.));
  calc_d_();
}

double Table1D::linear_interpolation(const double value0) const {
  const int n0 = num();
  TRACE("num0 " << n0);
  const int i0 = value0*(n0 - 1);
  TRACE("value0 " << value0);
  TRACE("i0 " << i0);
  int i02 = i0 + 1;
  if (i02 == n0) i02 = i0;
  TRACE("i02 " << i02);
  const double d0 = bin_spacing_;
  const double v0 = i0 * d0, vv0 = v0 + d0;
  TRACE("v0 " << v0);
  TRACE("vv0 " << vv0);
  const double xd0 = (value0 - v0) / (vv0 - v0);
  TRACE("xd0 " << xd0);
  TRACE("size0 " << data_.size());
  const double c00 = data_[i0] * (1-xd0) + xd0*data_[i02];
  TRACE("c000 " << data_[i0]);
  TRACE("c100 " << data_[i02]);
  TRACE("c00 " << c00);
  return c00;
}

double Table1D::forward_difference_interpolation(const double value0) const {
  const double sds = value0/bin_spacing_;
  TRACE("sds " << sds);
  const int k = int(sds);
  ASSERT(k + 2 < num(), "k: " << k << " beyond num: " << num());
  TRACE("k " << k);
  const double xi = sds - k;
  TRACE("xi " << xi);
  const double vk = data_[k];
  TRACE("vk " << vk);
  const double vk1 = data_[k + 1];
  TRACE("vk1 " << vk1);
  const double vk2 = data_[k + 2];
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

void Table1D::add(const Table1D& table) { feasst::add(table.data_, &data_); }

void Table2D::calc_d_() {
  bin_spacing_ = std::vector<double>({
    bin_spacing(num0()),
    bin_spacing(num1())});
}

Table2D::Table2D(argtype args) : Table2D(&args) { check_all_used(args); }
Table2D::Table2D(argtype * args) : Table() {
  const int num0 = integer("num0", args, 1);
  const int num1 = integer("num1", args, 1);
  resize(num0, num1, &data_);
  fill(dble("default_value", args, 0.), &data_);
  calc_d_();
}

double Table2D::linear_interpolation(const double value0,
    const double value1) const {
  const int n0 = num0();
  const int n1 = num1();
  TRACE("num0 " << n0);
  TRACE("num1 " << n1);
  const int i0 = value0*(n0 - 1);
  const int i1 = value1*(n1 - 1);
  TRACE("value0 " << value0);
  TRACE("value1 " << value1);
  TRACE("i0 " << i0 << " i1 " << i1);
  int i02 = i0 + 1, i12 = i1 + 1;
  if (i02 == n0) i02 = i0;
  if (i12 == n1) i12 = i1;
  TRACE("i02 " << i02 << " i12 " << i12);
  const double d0 = bin_spacing_[0];
  const double d1 = bin_spacing_[1];
  const double v0 = i0 * d0, vv0 = v0 + d0;
  const double v1 = i1 * d1, vv1 = v1 + d1;
  TRACE("v0 " << v0 << " v1 " << v1);
  TRACE("vv0 " << vv0 << " vv1 " << vv1);
  const double xd0 = (value0 - v0) / (vv0 - v0);
  const double xd1 = (value1 - v1) / (vv1 - v1);
  TRACE("xd0 " << xd0 << " xd1 " << xd1);
  TRACE("size0 " << data_.size());
  TRACE("size1 " << data_[0].size());
  const double c00 = data_[i0][i1] *(1-xd0) + xd0*data_[i02][i1];
  const double c10 = data_[i0][i12]*(1-xd0) + xd0*data_[i02][i12];
  TRACE("c000 " << data_[i0][i1] << " c010 " << data_[i0][i12]);
  TRACE("c100 " << data_[i02][i1] << " c110 " << data_[i02][i12]);
  TRACE("c00 " << c00 << " c10 " << c10);
  const double c0 = c00 * (1-xd1) + xd1*c10;
  TRACE("c0 " << c0);
  return c0;
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

void Table2D::add(const Table2D& table) { feasst::add(table.data_, &data_); }

void Table3D::calc_d_() {
  bin_spacing_ = std::vector<double>({
    bin_spacing(num0()),
    bin_spacing(num1()),
    bin_spacing(num2())});
}

Table3D::Table3D(argtype args) : Table3D(&args) { check_all_used(args); }
Table3D::Table3D(argtype * args) : Table() {
  const int num0 = integer("num0", args, 1);
  const int num1 = integer("num1", args, 1);
  const int num2 = integer("num2", args, 1);
  resize(num0, num1, num2, &data_);
  fill(dble("default_value", args, 0.), &data_);
  calc_d_();
}

double Table3D::linear_interpolation(const double value0,
    const double value1, const double value2) const {
  const int n0 = num0();
  TRACE("num0 " << n0);
  const int n1 = num1();
  TRACE("num1 " << n1);
  const int n2 = num2();
  TRACE("num2 " << n2);
  const int i0 = value0*(n0 - 1);
  const int i1 = value1*(n1 - 1);
  const int i2 = value2*(n2 - 1);
  TRACE("value0 " << value0);
  TRACE("value1 " << value1);
  TRACE("value2 " << value2);
  ASSERT(value0 >= 0 && value0 <= 1, "err");
  ASSERT(value1 >= 0 && value1 <= 1, "err");
  ASSERT(value2 >= 0 && value2 <= 1, "err");
  TRACE("i0 " << i0 << " i1 " << i1 << " i2 " << i2);
  int i02 = i0 + 1, i12 = i1 + 1, i22 = i2 + 1;
  if (i02 == n0) i02 = i0;
  if (i12 == n1) i12 = i1;
  if (i22 == n2) i22 = i2;
  TRACE("i02 " << i02 << " i12 " << i12 << " i22 " << i22);
  const double d0 = bin_spacing_[0];
  const double d1 = bin_spacing_[1];
  const double d2 = bin_spacing_[2];
  TRACE("d0 " << d0 << " d1 " << d1 << " d2 " << d2);
  const double v0 = i0 * d0, vv0 = v0 + d0;
  const double v1 = i1 * d1, vv1 = v1 + d1;
  const double v2 = i2 * d2, vv2 = v2 + d2;
  TRACE("v0 " << v0 << " v1 " << v1 << " v2 " << v2);
  TRACE("vv0 " << vv0 << " vv1 " << vv1 << " vv2 " << vv2);
  const double xd0 = (value0 - v0) / (vv0 - v0);
  const double xd1 = (value1 - v1) / (vv1 - v1);
  const double dv2 = vv2 - v2;
  double xd2;
  if (std::abs(dv2) < NEAR_ZERO) {
    xd2 = 0;
  } else {
    xd2 = (value2 - v2) / (vv2 - v2);
  }
  ASSERT(xd0 >= 0 && xd0 <= 1, "xd0 " << xd0);
  ASSERT(xd1 >= 0 && xd1 <= 1, "xd1 " << xd1);
  ASSERT(xd2 >= 0 && xd2 <= 1, "xd2 " << xd2);
  TRACE("xd0 " << xd0 << " xd1 " << xd1 << " xd2 " << xd2);
  TRACE("size0 " << data_.size());
  TRACE("size1 " << data_[0].size());
  TRACE("size2 " << data_[0][0].size());
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

void Table3D::serialize(std::ostream& ostr) const {
  feasst_serialize_version(6867, ostr);
  feasst_serialize(data_, ostr);
}

Table3D::Table3D(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6867, "version: " << version);
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

}  // namespace feasst
