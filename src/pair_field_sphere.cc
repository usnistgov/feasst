/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./pair_field_sphere.h"

namespace feasst {

PairFieldSphere::PairFieldSphere(Space * space)
  : Pair(space),
    PairField(space) {
  defaultConstruction_();
}

PairFieldSphere::PairFieldSphere(Space * space, const char* fileName)
  : Pair(space, fileName),
    PairField(space, fileName) {
  defaultConstruction_();
  radius_ = fstoi("radius", fileName);
  xSphere_ = fstoi("xSphere", fileName);
  ySphere_ = fstoi("ySphere", fileName);
  zSphere_ = fstoi("zSphere", fileName);
}

void PairFieldSphere::defaultConstruction_() {
  className_.assign("PairFieldSphere");
  initSphere();
  ASSERT(dimen_ > 1, "Assumes dimension > 1");
}

void PairFieldSphere::writeRestart(const char* fileName) {
  PairField::writeRestart(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  file << "# radius " << radius_ << endl;
  file << "# xSphere " << xSphere_ << endl;
  file << "# ySphere " << ySphere_ << endl;
  file << "# zSphere " << zSphere_ << endl;
}

void PairFieldSphere::fieldSite_(const int &siteType, double * energy, double * force,
  const double &x, const double &y, const double &z) {
  * energy = 0;
  const double dx = x - xSphere_,
               dy = y - ySphere_,
               dz = z - ySphere_,
               radial_dist = sqrt(dx*dx + dy*dy + dz*dz);
  fieldSiteSphere_(siteType, energy, force, radius_ - radial_dist);
}

void PairFieldSphere::initSphere(const double radius, const double xSphere,
  const double ySphere, const double zSphere) {
  radius_ = radius;
  xSphere_ = xSphere;
  ySphere_ = ySphere;
  zSphere_ = zSphere;
}

shared_ptr<PairFieldSphere> makePairFieldSphere(std::shared_ptr<Space> space) {
  return make_shared<PairFieldSphere>(space.get());
}

}  // namespace feasst
