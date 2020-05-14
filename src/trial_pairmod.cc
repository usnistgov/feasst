/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./trial_pairmod.h"
#include "./mc_wltmmc.h"

namespace feasst {

TrialPairMod::TrialPairMod() : Trial() {
  defaultConstruction_();
}
TrialPairMod::TrialPairMod(
  Pair *pair,
  Criteria *criteria)
  : Trial(pair, criteria) {
  defaultConstruction_();
}
TrialPairMod::TrialPairMod(const char* fileName,
  Pair *pair,
  Criteria *criteria)
  : Trial(pair, criteria, fileName) {
  defaultConstruction_();

  // initialize lnzRamp
  string strtmp = fstos("lnzMin", fileName);
  if (!strtmp.empty()) {
    initActivRamp(std::stod(strtmp), fstod("lnzMax", fileName));
  }
}

/**
 * write restart file
 */
void TrialPairMod::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  if (lnzRampOn_ == 1) {
    file << "# lnzMin " << lnzMin_ << endl;
    file << "# lnzMax " << lnzMax_ << endl;
  }
}

/**
 * default construction
 */
void TrialPairMod::defaultConstruction_() {
  className_.assign("TrialPairMod");
  trialType_.assign("move");
  verbose_ = 0;
  lnzRampOn_ = 0.;
}

/**
 * Attempt trial
 */
void TrialPairMod::attempt1_() {
  if (verbose_==1) {
    cout << std::setprecision(std::numeric_limits<double>::digits10+2)
         << "attempting gca " << pair_->peTot() << endl;
  }
  if (space()->nMol() <= 0) {
    trialMoveDecide_(0, 0);   // ensured rejection, however, criteria can update
    return void();
  }

  // store old state and energy
  trialMoveRecordAll_(0);
  const double orderOld = pair_->order();

  // attempt to vary the potential function
  //const double orderNew = orderOld + (uniformRanNum()-0.5)*maxMoveParam;
  double orderNew;
  if (uniformRanNum() < 0.5) {
    orderNew = orderOld + maxMoveParam;
  } else {
    orderNew = orderOld - maxMoveParam;
  }
  pair_->setOrder(orderNew);

  // if orderNew goes beyond certain bounds, it may be modified by the pair
  // class
  const double orderNewActual = pair_->order();

  // reject if bounds were hit
  double lnznew = -1;
  if (fabs(orderOld - orderNewActual) < DTOL) {
    reject_ = 1;
    de_ = 0;
    lnpMet_ = std::numeric_limits<double>::min();
  } else {
    reject_ = 0;
    de_ = pair_->allPartEnerForce(1) - peOld_;

    lnpMet_ = -criteria_->beta()*de_;
    // obtain new chemical potential
    if (lnzRampOn_) {
      ASSERT(0, "untested lnzRamp code");
      const double f = (orderNewActual - pair_->orderMin())
                       / (pair_->orderMax() - pair_->orderMin());
      lnznew = (lnzMax_ - lnzMin_)*f+lnzMin_;
      const double dlnz = lnznew - log( criteria_->activ(0) );
      lnpMet_ += space()->nMol()*dlnz;
    }

    // augment lnpMet by conjugate to shape
    if (orderNew > orderOld) {
      lnpMet_ += log(criteria_->activ(0));
    } else {
      lnpMet_ -= log(criteria_->activ(0));
    }
  }

  // WARNING: MAJOR HACK for bad interface
  // send the new order to criteria as a string in trialType_
  std::stringstream ss;
  ss << orderNewActual;

  if (criteria_->accept(lnpMet_, pair_->peTot() + de_, ss.str().c_str(),
                        reject_) == 1) {
    WARN(verbose_ == 1, "pairMod accepted " << de_);
    trialAccept_();
    pair_->update(space()->listAtoms(), 0, "update");
    if (lnzRampOn_) criteria_->activset( exp( lnznew ) );
    //cout << orderNewActual << " " << criteria_->activ()<< endl;
  } else {
    WARN(verbose_ == 1, "pairMod rejected " << de_);
    pair_->setOrder(orderOld);
    if (pair_->orderName().compare("MMsigma3") != 0) {
      space()->restoreAll();
      if (space()->cellType() > 0) space()->updateCellofallMol();
    }
    trialReject_();
  }
}

std::string TrialPairMod::printStat(const bool header) {
  std::stringstream stat;
  stat << Trial::printStat(header);
  if (header) {
    stat << "pairOrder ";
  } else {
    stat << pair_->order() << " ";
  }
  return stat.str();
}

shared_ptr<TrialPairMod> makeTrialPairMod(Pair *pair,
  Criteria *criteria) {
  return make_shared<TrialPairMod>(pair, criteria);
}

shared_ptr<TrialPairMod> makeTrialPairMod() {
  return make_shared<TrialPairMod>();
}

void pairModTrial(MC *mc, const double maxMoveParam) {
  shared_ptr<TrialPairMod> trial = make_shared<TrialPairMod>();
  if (maxMoveParam != -1) trial->maxMoveParam = maxMoveParam;
  mc->initTrial(trial);
}

void pairModTrial(shared_ptr<MC> mc, const double maxMoveParam) {
  pairModTrial(mc.get(), maxMoveParam);
}

void pairModTrial(WLTMMC *mc, double maxMoveParam) {
  shared_ptr<TrialPairMod> trial = make_shared<TrialPairMod>();
  // by default, the max move parameter is the size of the WL bin
  if (maxMoveParam == -1) {
    maxMoveParam = mc->c()->mBin();
    trial->maxMoveFlag = -1;
  }
  trial->maxMoveParam = maxMoveParam;
  mc->initTrial(trial);
}

void pairModTrial(shared_ptr<WLTMMC> mc, const double maxMoveParam) {
  pairModTrial(mc.get(), maxMoveParam);
}

}  // namespace feasst

