/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./pair_field_cylinder_sw.h"

namespace feasst {

PairFieldCylinderSW::PairFieldCylinderSW(Space * space)
  : Pair(space),
    PairField(space),
    PairFieldCylinder(space),
    PairFieldSW(space) {
  defaultConstruction_();
}

PairFieldCylinderSW::PairFieldCylinderSW(Space * space, const char* fileName)
  : Pair(space, fileName),
    PairField(space, fileName),
    PairFieldCylinder(space, fileName),
    PairFieldSW(space, fileName) {
  defaultConstruction_();
}

void PairFieldCylinderSW::defaultConstruction_() {
  className_.assign("PairFieldCylinderSW");
}

void PairFieldCylinderSW::writeRestart(const char* fileName) {
  PairFieldCylinder::writeRestart(fileName);
  PairFieldSW::writeRestart(fileName);
}

shared_ptr<PairFieldCylinderSW> makePairFieldCylinderSW(std::shared_ptr<Space> space) {
  return make_shared<PairFieldCylinderSW>(space.get());
}

}  // namespace feasst
