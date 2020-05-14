/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./pair_field_slit_sw.h"

namespace feasst {

PairFieldSlitSW::PairFieldSlitSW(Space * space)
  : Pair(space),
    PairField(space),
    PairFieldSlit(space),
    PairFieldSW(space) {
  defaultConstruction_();
}

PairFieldSlitSW::PairFieldSlitSW(Space * space, const char* fileName)
  : Pair(space, fileName),
    PairField(space, fileName),
    PairFieldSlit(space, fileName),
    PairFieldSW(space, fileName) {
  defaultConstruction_();
}

void PairFieldSlitSW::defaultConstruction_() {
  className_.assign("PairFieldSlitSW");
}

void PairFieldSlitSW::writeRestart(const char* fileName) {
  PairFieldSlit::writeRestart(fileName);
  PairFieldSW::writeRestart(fileName);
//  std::ofstream file(fileName, std::ios_base::app);
}

shared_ptr<PairFieldSlitSW> makePairFieldSlitSW(std::shared_ptr<Space> space) {
  return make_shared<PairFieldSlitSW>(space.get());
}

shared_ptr<PairFieldSlitSW> makePairFieldSlitSW(Space * space) {
  return make_shared<PairFieldSlitSW>(space);
}

}  // namespace feasst
