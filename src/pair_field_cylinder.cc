/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./pair_field_cylinder.h"

namespace feasst {

PairFieldCylinder::PairFieldCylinder(Space * space)
  : Pair(space),
    PairField(space) {
  defaultConstruction_();
}

PairFieldCylinder::PairFieldCylinder(Space * space, const char* fileName)
  : Pair(space, fileName),
    PairField(space, fileName) {
  defaultConstruction_();
  confine_dimen_ = fstoi("confineDimen", fileName);
  radius_ = fstoi("radius", fileName);
  xCylinder_ = fstoi("xCylinder", fileName);
  yCylinder_ = fstoi("yCylinder", fileName);
  zCylinder_ = fstoi("zCylinder", fileName);
}

void PairFieldCylinder::defaultConstruction_() {
  className_.assign("PairFieldCylinder");
  initCylinder();
  ASSERT(space()->dimen() == 3, "assumes 3D");
}

void PairFieldCylinder::writeRestart(const char* fileName) {
  PairField::writeRestart(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  file << "# confineDimen " << confine_dimen_ << endl;
  file << "# radius " << radius_ << endl;
  file << "# xCylinder " << xCylinder_ << endl;
  file << "# yCylinder " << yCylinder_ << endl;
  file << "# zCylinder " << zCylinder_ << endl;
}

void PairFieldCylinder::fieldSite_(const int &siteType, double * energy, double * force,
  const double &x, const double &y, const double &z) {
  * energy = 0;

  // find distance between point and cylindrical axis
  // http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
  // x1 and x2 are points along the axis, x0 is a test point
  vector<double> x0(3, 0), x1 = x0, x2 = x0;
  x0[0] = x; x0[1] = y; x0[2] = z;
  x1[0] = xCylinder_; x1[1] = yCylinder_; x1[2] = zCylinder_;
  x2 = x1;
  x2[confine_dimen_] += 1;
  const double radial_dist = closestDistancePointToAxis(x0, x1, x2);
  fieldSiteCylinder_(siteType, energy, force, radius_ - radial_dist);
}

void PairFieldCylinder::initCylinder(const int dimen, const double radius,
  const double xCylinder, const double yCylinder, const double zCylinder) {
  confine_dimen_ = dimen;
  radius_ = radius;
  xCylinder_ = xCylinder;
  yCylinder_ = yCylinder;
  zCylinder_ = zCylinder;
}

shared_ptr<PairFieldCylinder> makePairFieldCylinder(std::shared_ptr<Space> space) {
  return make_shared<PairFieldCylinder>(space.get());
}

}  // namespace feasst
