
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "utils/include/serialize_extra.h"
#include "math/include/recursive_table.h"

namespace feasst {

RecursiveTable1D::RecursiveTable1D(argtype * args) : Table1D(args) {
  nested_.resize(num() - 1);
}
RecursiveTable1D::RecursiveTable1D(argtype args) : RecursiveTable1D(&args) { feasst_check_all_used(args); }
RecursiveTable1D::~RecursiveTable1D() {}

void RecursiveTable1D::serialize(std::ostream& ostr) const {
  Table1D::serialize(ostr);
  feasst_serialize_version(1864, ostr);
  feasst_serialize(nested_, ostr);
}

RecursiveTable1D::RecursiveTable1D(std::istream& istr) : Table1D(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 1864, "version: " << version);
  feasst_deserialize(&nested_, istr);
}

void RecursiveTable1D::insert(const int bin, const RecursiveTable1D& nested) {
  ASSERT(bin >= 0 && bin < num(), "bin:" << bin << " must be >0 and < num:"
    << num());
  std::stringstream ss;
  nested.serialize(ss);
  nested_[bin] = std::make_unique<RecursiveTable1D>(ss);
}

double RecursiveTable1D::linear_interpolation(const double value0) const {
  const int n0 = num();
  TRACE("num0 " << n0);
  const int i0 = value0*(n0 - 1);
  TRACE("value0 " << value0);
  TRACE("i0 " << i0);
  int i02 = i0 + 1;
  if (i02 == n0) i02 = i0;
  TRACE("i02 " << i02);
  const double d0 = bin_spacing();
  const double v0 = i0 * d0, vv0 = v0 + d0;
  TRACE("v0 " << v0);
  TRACE("vv0 " << vv0);
  const double xd0 = (value0 - v0) / d0;
  TRACE("xd0 " << xd0);
  TRACE("size0 " << data().size());
  RecursiveTable1D * nested = nested_[i0].get();
  if (nested) {
    TRACE("begin nested");
    return nested->linear_interpolation(xd0);
  }
  const double c00 = data(i0) * (1-xd0) + xd0*data(i02);
  TRACE("c000 " << data(i0));
  TRACE("c100 " << data(i02));
  TRACE("c00 " << c00);
  return c00;
}

double RecursiveTable1D::forward_difference_interpolation(const double value0) const {
  FATAL("not implemented.");
}

}  // namespace feasst
