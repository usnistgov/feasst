/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./trial_swap.h"
#include "./mc.h"

namespace feasst {

TrialSwap::TrialSwap(const char* molTypeA,
  const char* molTypeB)
  : Trial(),
    molTypeA_(molTypeA),
    molTypeB_(molTypeB) {
  defaultConstruction_();
}

TrialSwap::TrialSwap(
  Pair *pair,
  Criteria *criteria,
  const char* molTypeA,
  const char* molTypeB)
  : Trial(pair, criteria),
    molTypeA_(molTypeA),
    molTypeB_(molTypeB) {
  defaultConstruction_();
}

TrialSwap::TrialSwap(const char* fileName,
  Pair *pair,
  Criteria *criteria)
  : Trial(pair, criteria, fileName) {
  defaultConstruction_();
  molTypeA_ = fstos("molTypeA", fileName);
  molTypeB_ = fstos("molTypeB", fileName);
}

void TrialSwap::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  file << "# molTypeA " << molTypeA_ << endl;
  file << "# molTypeB " << molTypeB_ << endl;
}

void TrialSwap::defaultConstruction_() {
  className_.assign("TrialSwap");
  trialType_.assign("move");
  verbose_ = 0;
}

void TrialSwap::attempt1_() {
  WARN(verbose_ == 1, "attempting to " << trialType_);
  const int nMolTypes = space()->addMolList().size();
  ASSERT(nMolTypes > 1,
    "xswap move requires nMolTypes(" << nMolTypes << ") > 1");

  // randomly choose a molecule
  int iMolOld = -1, iMolIndexOld = -1;
  string molTypeOld;
  if (space()->nMol() == 0) {
    reject_ = 1;
  } else {
    mpart_ = space()->randMol();
    iMolOld = space()->mol()[mpart_[0]];
    iMolIndexOld = space()->molid()[iMolOld];
    molTypeOld = space()->moltype()[iMolOld];
  }

  // check if it is of type A or B
  string molTypeNew;
  if (molTypeOld == molTypeA_) {
    molTypeNew = molTypeB_;
  } else if (molTypeOld == molTypeB_) {
    molTypeNew = molTypeA_;
  } else {
    reject_ = 1;
  }

  // record position and numbers of molecules
  vector<double> xAdd;
  int iMolIndexNew = -1;
  if (reject_ != 1) {
    iMolIndexNew = space()->findAddMolListIndex(molTypeNew);
    const int iMol = space()->randMolofType(iMolIndexOld);
    mpart_ = space()->imol2mpart(iMol);
    const int ipart = space()->mol2part()[iMol];
    for (int dim = 0; dim < space()->dimen(); ++dim) {
      xAdd.push_back(space()->x(ipart, dim));
    }
  }

  // record energy of molecule then swap
  if (reject_ != 1) {
    de_ = -1. * pair_->multiPartEner(mpart_, 0);

    // remove molecule
    pair_->delPart(mpart_);
    space()->delPart(mpart_);

    // add new molecule at the same position
    pair_->addMol(xAdd, molTypeNew.c_str());

    // record energy of new molecule
    mpart_ = space()->lastMolIDVec();
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
    trialAccept_();
    pair_->update(de_);

  // if not accepted, swap again
  } else {
    if (reject_ != 1) {
      pair_->delPart(mpart_);
      space()->delPart(mpart_);
      pair_->addMol(xAdd, molTypeOld.c_str());
    }
    trialReject_();
  }

  // record statistics on the composition
  if (reject_ != 1) {
    const int iMolIndexOld = space()->findAddMolListIndex(molTypeOld);
    const int nMolOld = space()->nMolType()[iMolIndexOld];
    const int iMolIndexNew = space()->findAddMolListIndex(molTypeNew);
    const int nMolNew = space()->nMolType()[iMolIndexNew];
    if (molTypeOld == molTypeA_) {
      nA_.accumulate(nMolOld);
      nB_.accumulate(nMolNew);
    } else if (molTypeOld == molTypeB_) {
      nA_.accumulate(nMolNew);
      nB_.accumulate(nMolOld);
    }
  }
}

string TrialSwap::printStat(const bool header) {
  stringstream stat;
  stat << Trial::printStat(header);
  if (header) {
    stat << "Na Nb ";
  } else {
    int iMolIndex = space()->findAddMolListIndex(molTypeA_);
    stat << space()->nMolType()[iMolIndex] << " ";
    iMolIndex = space()->findAddMolListIndex(molTypeB_);
    stat << space()->nMolType()[iMolIndex] << " ";
  }
  return stat.str();
}

shared_ptr<TrialSwap> makeTrialSwap(Pair *pair,
  Criteria *criteria, const char* molTypeA, const char* molTypeB) {
  return make_shared<TrialSwap>(pair, criteria, molTypeA, molTypeB);
}

shared_ptr<TrialSwap> makeTrialSwap(const char* molTypeA, const char* molTypeB) {
  return make_shared<TrialSwap>(molTypeA, molTypeB);
}

void swapTrial(MC *mc, const char* molTypeA, const char* molTypeB) {
  shared_ptr<TrialSwap> trial = make_shared<TrialSwap>(molTypeA, molTypeB);
  mc->initTrial(trial);
}

void swapTrial(shared_ptr<MC> mc, const char* molTypeA, const char* molTypeB) {
  swapTrial(mc.get(), molTypeA, molTypeB);
}

}  // namespace feasst

