#include "./trial_swap.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

TrialSwap::TrialSwap(const char* molTypeA,
  const char* molTypeB)
  : Trial(),
    molTypeA_(molTypeA),
    molTypeB_(molTypeB) {
  defaultConstruction();
}
TrialSwap::TrialSwap(Space *space,
  Pair *pair,
  Criteria *criteria,
  const char* molTypeA,
  const char* molTypeB)
  : Trial(space, pair, criteria),
    molTypeA_(molTypeA),
    molTypeB_(molTypeB) {
  defaultConstruction();
}
TrialSwap::TrialSwap(const char* fileName,
  Space *space,
  Pair *pair,
  Criteria *criteria)
  : Trial(space, pair, criteria, fileName) {
  defaultConstruction();
  molTypeA_ = fstos("molTypeA", fileName);
  molTypeB_ = fstos("molTypeB", fileName);
}

/**
 * write restart file
 */
void TrialSwap::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  file << "# molTypeA " << molTypeA_ << endl;
  file << "# molTypeB " << molTypeB_ << endl;
}

/**
 * default constructor
 */
void TrialSwap::defaultConstruction() {
  className_.assign("TrialSwap");
  trialType_.assign("move");
  verbose_ = 0;
}

/**
 * clone design pattern
 */
TrialSwap* TrialSwap::clone(Space* space, Pair *pair, Criteria *criteria) const {
  TrialSwap* t = new TrialSwap(*this);
  t->reconstruct(space, pair, criteria);
  return t;
}
shared_ptr<TrialSwap> TrialSwap::cloneShrPtr(
  Space* space, Pair* pair, Criteria* criteria) const {
  return(std::static_pointer_cast<TrialSwap, Trial>
    (cloneImpl(space, pair, criteria)));
}
shared_ptr<Trial> TrialSwap::cloneImpl(
  Space* space, Pair *pair, Criteria *criteria) const {
  shared_ptr<TrialSwap> t = make_shared<TrialSwap>(*this);
  t->reconstruct(space, pair, criteria);
  return t;
}


/**
 * Attempt to randomly swap molecule types
 */
void TrialSwap::attempt1() {
  WARN(verbose_ == 1, "attempting to " << trialType_);
  const int nMolTypes = space_->addMolList().size();
  ASSERT(nMolTypes > 1,
    "xswap move requires nMolTypes(" << nMolTypes << ") > 1");

  // randomly choose a molecule
  int iMolOld = -1, iMolIndexOld = -1;
  string molTypeOld;
  if (space_->nMol() == 0) {
    reject_ = 1;
  } else {
    mpart_ = space_->randMol();
    iMolOld = space_->mol()[mpart_[0]];
    iMolIndexOld = space_->molid()[iMolOld];
     molTypeOld = space_->moltype()[iMolOld];
  }

  //cout << "iMolOld " << iMolOld << endl;

  // check if it is of type A or B
  string molTypeNew;
  if (molTypeOld == molTypeA_) {
    molTypeNew = molTypeB_;
  } else if (molTypeOld == molTypeB_) {
    molTypeNew = molTypeA_;
  } else {
    reject_ = 1;
  }

  //cout << "molTypeOld " << molTypeOld << endl;
  //cout << "molTypeNew " << molTypeNew << endl;

  // record position and numbers of molecules
  vector<double> xAdd;
  int iMolIndexNew = -1;
  if (reject_ != 1) {
    iMolIndexNew = space_->findAddMolListIndex(molTypeNew);
    const int iMol = space_->randMolofType(iMolIndexOld);
    mpart_ = space_->imol2mpart(iMol);
    const int ipart = space_->mol2part()[iMol];
    for (int dim = 0; dim < space_->dimen(); ++dim) {
      xAdd.push_back(space_->x(ipart, dim));
    }
  }

  // record energy of molecule then swap
  if (reject_ != 1) {
    de_ = -1. * pair_->multiPartEner(mpart_, 0);

    // remove molecule
    pair_->delPart(mpart_);
    space_->delPart(mpart_);

    // add new molecule at the same position
    space_->xAdd = xAdd;
    space_->addMol(molTypeNew.c_str());
    pair_->addPart();

    // record energy of new molecule
    mpart_ = space_->lastMolIDVec();
    de_ += pair_->multiPartEner(mpart_, 0);

    lnpMet_ += -criteria_->beta()*de_
      + log(criteria_->activ(iMolIndexNew))
      - log(criteria_->activ(iMolIndexOld));
  }

  if (reject_ == 1) {
    preFac_ = 0;
    de_ = 0;
    lnpMet_ = std::numeric_limits<double>::min();
  }

  if (verbose_ == 1) {
    cout << "de " << de_ << " rej " << reject_ << " lpm " << lnpMet_ <<  endl;
  }

  // acceptance criteria
  if (criteria_->accept(lnpMet_, pair_->peTot() + de_,
                        trialType_.c_str(), reject_) == 1) {
    trialAccept();
    pair_->update(de_);

  // if not accepted, swap again
  } else {
    if (reject_ != 1) {
      pair_->delPart(mpart_);
      space_->delPart(mpart_);
      space_->xAdd = xAdd;
      space_->addMol(molTypeOld.c_str());
      pair_->addPart();
    }
    trialReject();
  }

  // record statistics on the composition
  if (reject_ != 1) {
    const int iMolIndexOld = space_->findAddMolListIndex(molTypeOld);
    const int nMolOld = space_->nMolType()[iMolIndexOld];
    const int iMolIndexNew = space_->findAddMolListIndex(molTypeNew);
    const int nMolNew = space_->nMolType()[iMolIndexNew];
    if (molTypeOld == molTypeA_) {
      nA_.accumulate(nMolOld);
      nB_.accumulate(nMolNew);
    } else if (molTypeOld == molTypeB_) {
      nA_.accumulate(nMolNew);
      nB_.accumulate(nMolOld);
    }
  }
}

/*
 * return string for status of trial
 */
string TrialSwap::printStat(const bool header) {
  stringstream stat;
  stat << Trial::printStat(header);
  if (header) {
    stat << "Na Nb ";
  } else {
    int iMolIndex = space_->findAddMolListIndex(molTypeA_);
    stat << space_->nMolType()[iMolIndex] << " ";
    iMolIndex = space_->findAddMolListIndex(molTypeB_);
    stat << space_->nMolType()[iMolIndex] << " ";
  }
  return stat.str();
}

///**
// * Attempt to randomly swap molecule types
// */
//void TrialSwap::attempt2() {
//  WARN(verbose_ == 1, "attempting to " << trialType_);
//  const int nMolTypes = space_->addMolList().size();
//  ASSERT(nMolTypes > 1,
//    "xswap move requires nMolTypes(" << nMolTypes << ") > 1");
//
//  // randomly choose to select an A or a B molecule
//  string molTypeOld = molTypeA_;
//  string molTypeNew = molTypeB_;
//  if (uniformRanNum() > 0.5) {
//    molTypeOld = molTypeB_;
//    molTypeNew = molTypeA_;
//  }
//  //cout << "ml " << molType << " " << molTypeNew << endl;
//  // choose a particle of this type and record its position
//  vector<double> xAdd;
//  const int iMolIndexOld = space_->findAddMolListIndex(molTypeOld);
//  const int nMolOld = space_->nMolType()[iMolIndexOld];
//  const int iMolIndexNew = space_->findAddMolListIndex(molTypeNew);
//  const int nMolNew = space_->nMolType()[iMolIndexNew];
//  if (nMolOld > 0) {
//    const int iMol = space_->randMolofType(iMolIndexOld);
//    mpart_ = space_->imol2mpart(iMol);
//    const int ipart = space_->mol2part()[iMol];
//    for (int dim = 0; dim < space_->dimen(); ++dim) {
//      xAdd.push_back(space_->x(ipart, dim));
//    }
//  } else {
//    reject_ = 1;
//  }
//
//  if (reject_ != 1) {
//    // record energy of molecule
//    de_ = -1. * pair_->multiPartEner(mpart_, 0);
//
//    // remove molecule
//    pair_->delPart(mpart_);
//    space_->delPart(mpart_);
//
//    // add new molecule at the same position
//    space_->xAdd = xAdd;
//    space_->addMol(molTypeNew.c_str());
//    pair_->addPart();
//
//    // record energy of new molecule
//    mpart_ = space_->lastMolIDVec();
//    de_ += pair_->multiPartEner(mpart_, 0);
//
//    lnpMet_ += -criteria_->beta()*de_
//      + log(criteria_->activ(iMolIndexNew))
//      - log(criteria_->activ(iMolIndexOld))
//      + log(static_cast<double>(nMolOld)/
//            static_cast<double>(nMolNew+1));
//  }
//
//  if (reject_ == 1) {
//    preFac_ = 0;
//    de_ = 0;
//    lnpMet_ = std::numeric_limits<double>::min();
//  }
//
//  if (verbose_ == 1) {
//    cout << "de " << de_ << " rej " << reject_ << " lpm " << lnpMet_ <<  endl;
//  }
//
//  // acceptance criteria
//  if (criteria_->accept(lnpMet_, pair_->peTot() + de_,
//                        trialType_.c_str(), reject_) == 1) {
//    trialAccept();
//    pair_->update(de_);
//
//  // if not accepted, swap again
//  } else {
//    if (reject_ != 1) {
//      pair_->delPart(mpart_);
//      space_->delPart(mpart_);
//      space_->xAdd = xAdd;
//      space_->addMol(molTypeOld.c_str());
//      pair_->addPart();
//    }
//    trialReject();
//  }
//}
//
//
//
//
//

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

