/**
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of this agreement (see LICENSE.txt) and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#include "./trial_xswap.h"
#include "./mc.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

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

/**
 * default construction
 */
void TrialXSwap::defaultConstruction_() {
  className_.assign("TrialXSwap");
  trialType_.assign("move");
  verbose_ = 0;
}

/**
 * Attempt trial
 */
void TrialXSwap::attempt1_() {
  if (verbose_ == 1) {
    cout << std::setprecision(std::numeric_limits<double>::digits10+2)
         << "attempting gca " << pair_->peTot() << endl;
  }
  if (space()->nMol() <= 1) {
    trialMoveDecide_(0, 0);   // ensured rejection, however, criteria can update
    return void();
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

  if ( (iMol == -1) || (jMol == -1) ) {
    reject_ = 1;
  } else {
    vector<int> ipart = space()->imol2mpart(iMol),
                jpart = space()->imol2mpart(jMol);
//        mpart_.reserve(ipart.size() + jpart.size());
//        mpart_.insert(mpart_.end(), ipart.begin(), ipart.end());
//        mpart_.insert(mpart_.end(), jpart.begin(), jpart.end());

    // peOld_ = pair_->multiPartEner(mpart_, 0);
    // pair_->update(mpart_, 0, "store");
    // space()->xStore(mpart_);
    // trialMoveRecord_();
    peOld_ = pair_->multiPartEner(ipart, 0) + pair_->multiPartEner(jpart, 0);
    space()->swapPositions(iMol, jMol);
    double peNew;
    peNew = pair_->multiPartEner(ipart, 0) + pair_->multiPartEner(jpart, 0);
    // cout << "peNew " << peNew << endl;
    // peNew = pair_->multiPartEner(mpart_, 0);
    // cout << "peNew " << peNew << endl;
    de_ = peNew - peOld_;
    lnpMet_ = -criteria_->beta()*de_;
    reject_ = 0;
  }

  // accept or reject
  if (criteria_->accept(lnpMet_, pair_->peTot() + de_, trialType_.c_str(),
                        reject_) == 1) {
    // pair_->update(mpart_, 0, "update");
    pair_->update(de_);
    trialAccept_();
  } else {
    // space()->restore(mpart_);
    if (reject_ != 1) space()->swapPositions(iMol, jMol);
    trialReject_();
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

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

