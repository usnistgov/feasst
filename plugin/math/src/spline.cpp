
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/spline.h"

namespace feasst {

LinearSpline1D::LinearSpline1D(const std::vector<double>& ycoord) {
  values_ = ycoord;
  dz_ = 1./static_cast<double>(values_.size() - 1);
}

void LinearSpline1D::serialize(std::ostream& sstr) const {
  feasst_serialize_version(6285, sstr);
  feasst_serialize(values_, sstr);
  feasst_serialize(dz_, sstr);
}

LinearSpline1D::LinearSpline1D(std::istream& sstr) {
  const int version = feasst_deserialize_version(sstr);
  ASSERT(version == 6285, "unrecognized verison: " << version);
  feasst_deserialize(&values_, sstr);
  feasst_deserialize(&dz_, sstr);
}

double LinearSpline1D::interpolate(const double z) const {
  ASSERT(z >= 0. && z <= 1., "z:" << z);
  if (z == 0) {
    return values_[0];
  } else if (z == 1) {
    return values_.back();
  }
  const int bin = int(z*static_cast<int>(values_.size() - 1));
  const double zz = z/dz_ - bin;
  const double c0 = values_[bin],
               c1 = values_[bin + 1];
  return c0 + zz*(c1 - c0);
}

LinearSpline2D::LinearSpline2D(const vec2& ycoord) {
  values_ = ycoord;
  dx_ = 1./static_cast<double>(values_.size() - 1);
  dy_ = 1./static_cast<double>(values_[0].size() - 1);
}

void LinearSpline2D::serialize(std::ostream& sstr) const {
  feasst_serialize_version(3240, sstr);
  feasst_serialize(values_, sstr);
  feasst_serialize(dx_, sstr);
  feasst_serialize(dy_, sstr);
}

LinearSpline2D::LinearSpline2D(std::istream& sstr) {
  const int version = feasst_deserialize_version(sstr);
  ASSERT(version == 3240, "unrecognized verison: " << version);
  feasst_deserialize(&values_, sstr);
  feasst_deserialize(&dx_, sstr);
  feasst_deserialize(&dy_, sstr);
}

double LinearSpline2D::interpolate(const double xs, const double ys) const {
  ASSERT(xs >= 0. && xs <= 1., "xs:" << xs);
  ASSERT(ys >= 0. && ys <= 1., "ys:" << ys);
//  if (xs == 0 || ys == 0) {
//    FATAL("imp");
//    //return values_[0];
//  } else if (xs == 1 || ys == 1) {
//    FATAL("imp");
//    //return values_.back();
//  }
  const int xbin = int(xs*static_cast<int>(values_.size() - 1));
  const int ybin = int(ys*static_cast<int>(values_[0].size() - 1));
  const double xx = xs/dx_ - xbin;
  const double yy = ys/dy_ - ybin;
  const double c00 = values_[xbin][ybin];
  const double c01 = values_[xbin][ybin + 1];
  const double c10 = values_[xbin + 1][ybin];
  const double c11 = values_[xbin + 1][ybin + 1];
  const double A = c00, B1 = c10 - A, B2 = c01 - A, B12 = c11 - c10 - B2;
  return A + B1*xx + B2*yy + B12*xx*yy;
}

QuadraticSpline1D::QuadraticSpline1D(const std::vector<double>& ycoord,
    const double right_endpoint_derivative) {
  values_ = ycoord;
  derivs_.resize(values_.size());
  derivs_.back() = right_endpoint_derivative;
  for (int i = static_cast<int>(derivs_.size() - 2); i >= 0; --i) {
    derivs_[i] = 2*(values_[i+1] - values_[i]) - derivs_[i + 1];
  }
  dz_ = 1./static_cast<double>(values_.size() - 1);
}

void QuadraticSpline1D::serialize(std::ostream& sstr) const {
  feasst_serialize_version(1674, sstr);
  feasst_serialize(values_, sstr);
  feasst_serialize(derivs_, sstr);
  feasst_serialize(dz_, sstr);
}

QuadraticSpline1D::QuadraticSpline1D(std::istream& sstr) {
  const int version = feasst_deserialize_version(sstr);
  ASSERT(version == 1674, "unrecognized verison: " << version);
  feasst_deserialize(&values_, sstr);
  feasst_deserialize(&derivs_, sstr);
  feasst_deserialize(&dz_, sstr);
}

double QuadraticSpline1D::interpolate(const double z) const {
  ASSERT(z >= 0. && z <= 1., "z:" << z);
  if (z == 0) {
    return values_[0];
  } else if (z == 1) {
    return values_.back();
  }
  const int size = static_cast<int>(derivs_.size());
  const int bin = int(z*(size - 1));
  const double dz = 1./static_cast<double>(size - 1);
  const double zz = z/dz - bin;
  const double lower = values_[bin];
  const double upper = values_[bin + 1];
  const double deriv = derivs_[bin + 1];
  const double dy = upper - lower;
  return lower + zz*((2*dy - deriv) + zz*(-dy + deriv));
}

}  // namespace feasst
