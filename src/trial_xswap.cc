/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./trial_xswap.h"
#include "./mc.h"

namespace feasst {

TrialXSwap::TrialXSwap() : Trial() {
  defaultConstruction_();
}

TrialXSwap::TrialXSwap(
  Pair *pair,
  Criteria *criteria)
  : Trial(pair, criteria) {
  defaultConstruction_();
}

TrialXSwap::TrialXSwap(const char* fileName,
  Pair *pair,
  Criteria *criteria)
  : Trial(pair, criteria, fileName) {
  defaultConstruction_();
}

void TrialXSwap::defaultConstruction_() {
  className_.assign("TrialXSwap");
  trialType_.assign("move");
  verbose_ = 0;
}

void TrialXSwap::attempt1_() {
  if (space()->nMol() <= 1) {
    trialMoveDecide_(0, 0);   // ensured rejection, however, criteria can update
    return;
  }

  // check that there is more than one type of molecule
  const int nMolTypes = space()->addMolList().size();
  ASSERT(nMolTypes > 1,
    "xswap move requires nMolTypes(" << nMolTypes << ") >= 1");

  // randomly pick two different types of molecules
  const int iType = uniformRanNum(0, nMolTypes-1);
  int iter = 0, jType, nmax = 1e6;
  bool term = false;
  while (!term && (iter < nmax)) {
    jType = uniformRanNum(0, nMolTypes-1);
    if (iType != jType) term = true;
    ++iter;
    if (iter == nmax) reject_ = 1;
  }

  // select a random molecule of each type
  const int iMol = space()->randMolofType(iType),
            jMol = space()->randMolofType(jType);

  WARN(verbose_ == 1, "swapping iMol " << iMol << " and jMol " << jMol);
  if ( (iMol == -1) || (jMol == -1) || (reject_ == 1) ) {
    reject_ = 1;
    trialMoveDecide_(0, 0);
    return;
  }

  mpart_ = space()->imol2mpart(iMol);
  const vector<int> mpart2 = space()->imol2mpart(jMol);
  mpart_.insert(mpart_.end(), mpart2.begin(), mpart2.end());
  std::sort(mpart_.begin(), mpart_.end());

  trialMoveRecord_();
  space()->swapPositions(iMol, jMol);
  if (space()->cellType() > 0) {
    space()->updateCellofiMol(iMol);
    space()->updateCellofiMol(jMol);
  }
  reject_ = 0;
  trialMoveDecide_(0, 1);
  if (space()->cellType() > 0) {
    space()->updateCellofiMol(iMol);
    space()->updateCellofiMol(jMol);
  }
}

shared_ptr<TrialXSwap> makeTrialXSwap(Pair *pair,
  Criteria *criteria) {
  return make_shared<TrialXSwap>(pair, criteria);
}

shared_ptr<TrialXSwap> makeTrialXSwap() {
  return make_shared<TrialXSwap>();
}

void xswapTrial(MC *mc) {
  shared_ptr<TrialXSwap> trial = make_shared<TrialXSwap>();
  mc->initTrial(trial);
}
void xswapTrial(shared_ptr<MC> mc) {
  xswapTrial(mc.get());
}

}  // namespace feasst

