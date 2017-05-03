#include "./trial_xswap.h"

namespace feasst {

TrialXSwap::TrialXSwap() : Trial() {
  defaultConstruction();
}
TrialXSwap::TrialXSwap(Space *space,
  Pair *pair,
  Criteria *criteria)
  : Trial(space, pair, criteria) {
  defaultConstruction();
}
TrialXSwap::TrialXSwap(const char* fileName,
  Space *space,
  Pair *pair,
  Criteria *criteria)
  : Trial(space, pair, criteria, fileName) {
  defaultConstruction();
}

/**
 * default construction
 */
void TrialXSwap::defaultConstruction() {
  className_.assign("TrialXSwap");
  trialType_.assign("move");
  verbose_ = 0;
}

/**
 * Attempt trial
 */
void TrialXSwap::attempt1() {
  if (verbose_ == 1) {
    cout << std::setprecision(std::numeric_limits<double>::digits10+2)
         << "attempting gca " << pair_->peTot() << endl;
  }
  if (space_->nMol() <= 1) {
    trialMoveDecide(0, 0);   // ensured rejection, however, criteria can update
    return void();
  }

  // check that there is more than one type of molecule
  const int nMolTypes = space_->addMolList().size();
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
  const int iMol = space_->randMolofType(iType),
            jMol = space_->randMolofType(jType);

  if ( (iMol == -1) || (jMol == -1) ) {
    reject_ = 1;
  } else {
    vector<int> ipart = space_->imol2mpart(iMol),
                jpart = space_->imol2mpart(jMol);
//        mpart_.reserve(ipart.size() + jpart.size());
//        mpart_.insert(mpart_.end(), ipart.begin(), ipart.end());
//        mpart_.insert(mpart_.end(), jpart.begin(), jpart.end());

    // peOld_ = pair_->multiPartEner(mpart_, 0);
    // pair_->update(mpart_, 0, "store");
    // space_->xStore(mpart_);
    // trialMoveRecord();
    peOld_ = pair_->multiPartEner(ipart, 0) + pair_->multiPartEner(jpart, 0);
    space_->swapPositions(iMol, jMol);
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
    trialAccept();
  } else {
    // space_->restore(mpart_);
    if (reject_ != 1) space_->swapPositions(iMol, jMol);
    trialReject();
  }
}

}  // namespace feasst

