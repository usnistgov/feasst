/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./trial_beta.h"
#include "./mc_wltmmc.h"

namespace feasst {

TrialBeta::TrialBeta()
  : Trial() {
  defaultConstruction_();
}
TrialBeta::TrialBeta(
  Pair *pair,
  Criteria *criteria)
  : Trial(pair, criteria) {
  defaultConstruction_();
}
TrialBeta::TrialBeta(const char* fileName,
  Pair *pair,
  Criteria *criteria)
  : Trial(pair, criteria, fileName) {
  defaultConstruction_();
}

/**
 * default construction
 */
void TrialBeta::defaultConstruction_() {
  className_.assign("TrialBeta");
  trialType_.assign("move");
  verbose_ = 0;
}

/**
 * Attempt trial
 */
void TrialBeta::attempt1_() {
  if (verbose_ ==1) {
    cout << std::setprecision(std::numeric_limits<double>::digits10+2)
         << "attempting gca " << pair_->peTot() << endl;
  }
  if (space()->nMol() <= 0) {
    trialMoveDecide_(0, 0);   // ensured rejection, however, criteria can update
    return void();
  }

  // randomly attempt to increase or decrease beta by maxMoveParam
  double dBeta = maxMoveParam;
  if (uniformRanNum() < 0.5) dBeta *= -1;
  lnpMet_ = -dBeta * pair_->peTot();

  // WARNING: MAJOR HACK for bad interface
  // send the new order to criteria as a string in trialType_
  std::stringstream ss;
  ss << criteria_->beta() + dBeta;

  if (criteria_->accept(lnpMet_, pair_->peTot(), ss.str().c_str(),
                        reject_) == 1) {
    WARN(verbose_ == 1, "beta transform accepted " << de_);
    trialAccept_();
    criteria_->betaset(criteria_->beta() + dBeta);
  } else {
    WARN(verbose_ == 1, "beta transform rejected " << de_);
    trialReject_();
  }
}

std::string TrialBeta::printStat(const bool header) {
  std::stringstream stat;
  stat << Trial::printStat(header);
  if (header) {
    stat << "beta ";
  } else {
    stat << criteria_->beta() << " ";
  }
  return stat.str();
}

shared_ptr<TrialBeta> makeTrialBeta(Pair *pair,
  Criteria *criteria) {
  return make_shared<TrialBeta>(pair, criteria);
}

shared_ptr<TrialBeta> makeTrialBeta() {
  return make_shared<TrialBeta>();
}

void betaTrial(MC *mc, const double maxMoveParam) {
  shared_ptr<TrialBeta> trial = make_shared<TrialBeta>();
  if (maxMoveParam != -1) trial->maxMoveParam = maxMoveParam;
  mc->initTrial(trial);
}
void betaTrial(shared_ptr<MC> mc, const double maxMoveParam) {
  betaTrial(mc.get(), maxMoveParam);
}
void betaTrial(WLTMMC *mc, double maxMoveParam) {
  shared_ptr<TrialBeta> trial = make_shared<TrialBeta>();
  // by default, the max move parameter is the size of the WL bin
  if (maxMoveParam == -1) {
    maxMoveParam = mc->c()->mBin();
  }
  trial->maxMoveParam = maxMoveParam;
  mc->initTrial(trial);
}
void betaTrial(shared_ptr<WLTMMC> mc, const double maxMoveParam) {
  betaTrial(mc.get(), maxMoveParam);
}

}  // namespace feasst

