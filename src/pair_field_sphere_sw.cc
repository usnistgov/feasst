/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./pair_field_sphere_sw.h"

namespace feasst {

PairFieldSphereSW::PairFieldSphereSW(Space * space)
  : Pair(space),
    PairField(space),
    PairFieldSphere(space),
    PairFieldSW(space) {
  defaultConstruction_();
}

PairFieldSphereSW::PairFieldSphereSW(Space * space,
  const char* fileName)
  : Pair(space, fileName),
    PairField(space, fileName),
    PairFieldSphere(space, fileName),
    PairFieldSW(space, fileName) {
  defaultConstruction_();
}

void PairFieldSphereSW::defaultConstruction_() {
  className_.assign("PairFieldSphereSW");
}

void PairFieldSphereSW::writeRestart(const char* fileName) {
  PairFieldSphere::writeRestart(fileName);
  PairFieldSW::writeRestart(fileName);
}

shared_ptr<PairFieldSphereSW> makePairFieldSphereSW(std::shared_ptr<Space> space) {
  return make_shared<PairFieldSphereSW>(space.get());
}

}  // namespace feasst
