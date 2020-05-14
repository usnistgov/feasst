/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./pair_field_slit_linear.h"

namespace feasst {

PairFieldSlitLinear::PairFieldSlitLinear(Space * space)
  : Pair(space),
    PairField(space),
    PairFieldSlit(space),
    PairFieldLinear(space) {
  defaultConstruction_();
}

PairFieldSlitLinear::PairFieldSlitLinear(Space * space, const char* fileName)
  : Pair(space, fileName),
    PairField(space, fileName),
    PairFieldSlit(space, fileName),
    PairFieldLinear(space, fileName) {
  defaultConstruction_();
}

void PairFieldSlitLinear::defaultConstruction_() {
  className_.assign("PairFieldSlitLinear");
}

void PairFieldSlitLinear::writeRestart(const char* fileName) {
  PairFieldSlit::writeRestart(fileName);
  PairFieldLinear::writeRestart(fileName);
//  std::ofstream file(fileName, std::ios_base::app);
}

shared_ptr<PairFieldSlitLinear> makePairFieldSlitLinear(std::shared_ptr<Space> space) {
  return make_shared<PairFieldSlitLinear>(space.get());
}

shared_ptr<PairFieldSlitLinear> makePairFieldSlitLinear(Space * space) {
  return make_shared<PairFieldSlitLinear>(space);
}

}  // namespace feasst
