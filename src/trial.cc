/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./trial.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

Trial::Trial() {
  defaultConstruction_();
}

Trial::Trial(Pair *pair, Criteria *criteria, const argtype &args)
  : pair_(pair),
    criteria_(criteria) {
  defaultConstruction_();
  argparse_.initArgs(className_, args);

  // parse maxMoveParam
  if (!argparse_.key("maxMoveParam").empty()) {
    maxMoveParam = stod(argparse_.str());
  }
}

Trial::Trial(Pair *pair, Criteria *criteria, const char* fileName)
  : pair_(pair),
    criteria_(criteria) {
  ASSERT(fileExists(fileName),
         "restart file(" << fileName << ") doesn't exist");
  defaultConstruction_();

  // initialize random number generator
  initRNG(fileName);

  accepted_ = fstoll("accepted", fileName);
  attempted_ = fstoll("attempted", fileName);
  string strtmp = fstos("numfirstbead", fileName);
  if (!strtmp.empty()) {
    nf_ = fstoi("numfirstbead", fileName);
    numFirstBeads(nf_);
  }
  const string rabovestr = fstos("rAbove", fileName);
  if (!rabovestr.empty()) {
    rAbove_ = stod(rabovestr);
    rBelow_ = stod(fstos("rBelow", fileName));
    initAVB(rAbove_, rBelow_);
  }

  strtmp = fstos("maxMoveParam", fileName);
  if (!strtmp.empty()) {
    maxMoveParam = stod(strtmp);
  }

  strtmp = fstos("confineUpper", fileName);
  if (!strtmp.empty()) {
    double lower = fstod("confineLower", fileName);
    int dim = fstoi("confineDim", fileName);
    confine(stod(strtmp), lower, dim);
  }

  strtmp = fstos("nPartTarget", fileName);
  if (!strtmp.empty()) {
    initializeNMolSeek(stoi(strtmp));
  }
}

void Trial::defaultConstruction_() {
  className_.assign("Trial");
  zeroStat();
  maxMoveParam = maxMoveParamDefault_;
  maxMoveFlag = 0;
  verbose_ = 0;
  nf_ = 0;
  rAbove_ = -1;
  rBelow_ = -1;
  avbOn_ = false;
  confineFlag_ = 0;
  initializeNMolSeek();
}

void Trial::reconstruct(Pair *pair, Criteria *criteria) {
  pair_ = pair;
  criteria_ = criteria;
  Base::reconstruct();
}

void Trial::zeroStat() {
  de_ = 0;
  deTot_ = 0;
  accepted_ = 0;
  attempted_ = 0;
}

void Trial::attempt() {
  ++attempted_;
  criteria_->store(pair_);
  preFac_ = 0;
  lnpMet_ = 0;
  reject_ = 0;
  de_ = 0.;
  def_ = 0.;
  attempt1_();
}

void Trial::trialAccept_() {
  ++accepted_;
  deTot_ += de_;
  WARN(verbose_ == 1, "accepted " << de_);
}

void Trial::trialReject_() {
  WARN(verbose_ == 1, "rejected " << de_);
  de_ = 0;
}

void Trial::numFirstBeads(const int nf) {
  ASSERT((nf == 0) || (nf > 1),
         "nf(" << nf << ") should be 0 or greater than 1");
  nf_ = nf;
  en_.resize(nf_);
  w_.resize(nf_);
  cpdf_.resize(nf_);
}

void Trial::initAVB(const double rAbove, const double rBelow) {
  ASSERT(rAbove > rBelow, "Null aggregation volume when rAbove = "
         << rAbove << " and rBelow = " << rBelow);
  ASSERT(rAbove <= 0.5*space()->minl(), "aggregation volume upper radius("
    << rAbove << ") cannot extend beyond periodic boundary conditions when"
    << "minimum box length is " << space()->minl());
  ASSERT(rAbove <= pair_->rCut(), "aggregation volume upper radius("
    << rAbove << ") cannot extend beyond pair cut-off " << pair_->rCut());
  rAbove_ = rAbove;
  rBelow_ = rBelow;
  vIn_ = volShell(rAbove_, rBelow_, space()->dimen());
  avbOn_ = true;
}

double Trial::multiFirstBead_(const int flag) {
  ASSERT((flag == 0) || (flag == 1) || (flag == 2) || (flag == 3),
    "Unrecognized flag for multiple first beads");

  // record position
  space()->xStoreMulti(mpart_, -1);
  pair_->cheapEnergy(1);
  for (int i = 0; i < nf_; ++i) {
    if (i != 0) {
      if (avbOn_) {
        space()->avb(mpart_, tmpart_, rAbove_, rBelow_, region_.c_str());
      } else {
        // -1 flag for displacement by half box length in any direction
        space()->randDisp(mpart_, -1);
        // -1 flag for completely random displacement
        space()->randRotate(mpart_, -1);
      }
      space()->wrap(mpart_);
      // -2 flag to store multiple positions
      space()->xStoreMulti(mpart_, -2);
    }
    en_[i] = pair_->multiPartEner(mpart_, flag);
    w_[i] = exp(-criteria_->beta()*en_[i]);
    // cout << "w" << i << " " << w_[i] << " " << en_[i] << " " << flag << endl;
    if (i != 0) w_[i] += w_[i-1];
  }
  const double wr = w_.back();
  // cout << "wr " << wr << endl;
  pair_->cheapEnergy(0);

  // select trial position from rosenblueth factor weights
  if (wr != 0) {
    // if new configuration, select f and multiply by weight
    if ( (flag == 1) || (flag == 3) ) {
      for (int i = 0; i < nf_; ++i) {
        cpdf_[i] = w_[i]/wr;
      }
      int f = ranFromCPDF(cpdf_);
      def_ = en_[f];
      space()->xStoreMulti(mpart_, f);
      return wr/static_cast<double>(nf_);

    // if old configuration, select first trial and divide by weight
    } else {
      def_ = en_[0];
      space()->xStoreMulti(mpart_, 0);
      // cout << "def_ " << def_ << " nf/wr " << nf_/wr << endl;
      return static_cast<double>(nf_)/wr;
    }
  } else {
    def_ = 0;
    return 0;
  }
}

void Trial::replaceCriteria(Criteria *criteria) {
  criteriaOld_ = criteria_;
  criteria_ = criteria;
}
void Trial::restoreCriteria() {
  if (criteriaOld_ != NULL) {
    criteria_ = criteriaOld_;
    criteriaOld_ = NULL;
  } else {
    ASSERT(0, "attempting to restore criteria, but not old criteria recorded");
  }
}

void Trial::writeRestartBase(const char* fileName) {
  fileBackUp(fileName);
  std::ofstream file(fileName);
  file << "# className " << className_ << endl;
  file << "# accepted " << accepted_ << endl;
  file << "# attempted " << attempted_ << endl;
  if (nf_ != 0) file << "# numfirstbead " << nf_ << endl;
  if (avbOn_) {
    file << "# rAbove " << rAbove_ << endl;
    file << "# rBelow " << rBelow_ << endl;
  }
  if (maxMoveParam != maxMoveParamDefault_) {
    file << "# maxMoveParam " << maxMoveParam << endl;
  }
  if (confineFlag_ == 1) {
    file << "# confineUpper " << confineUpper_ << endl;
    file << "# confineLower " << confineLower_ << endl;
    file << "# confineDim " << confineDim_ << endl;
  }
  if (nPartTarget_ != -1) {
    file << "# nPartTarget " << nPartTarget_ << endl;
  }

  // write random number generator state
  writeRngRestart(fileName);
}

double Trial::acceptPer() const {
  if (attempted_ != 0) {
    return static_cast<double>(accepted_)/attempted_;
  } else {
    return 0;
  }
}

void Trial::trialMoveRecord_() {
  if (criteria_->className() != "CriteriaMayer") {
    peOld_ = pair_->multiPartEner(mpart_, 0);
    pair_->update(mpart_, 0, "store");
  } else {
    peOld_ = pair_->peTot();
  }
  space()->xStore(mpart_);
}

void Trial::trialMoveRecordAll_(const int flag) {
  peOld_ = pair_->allPartEnerForce(flag);
  pair_->update(space()->listAtoms(), 0, "store");
  space()->xStoreAll();
}

void Trial::trialMoveDecide_(const double def,
  const double preFac) {
  if (preFac != 0) {
    // compute energy contribution of selected molecule in new configuration
    if (space()->cellType() > 0) {
      space()->updateCellofiMol(space()->mol()[mpart_.front()]);
    }
    const double pe = pair_->multiPartEner(mpart_, 1);
    de_ = pe - peOld_;
    lnpMet_ = log(preFac) - criteria_->beta()*(de_ - def);
    reject_ = 0;
    // if accepted, update
    if (criteria_->accept(lnpMet_, pair_->peTot() + de_, trialType_.c_str(),
                          reject_) == 1) {
      space()->wrap(mpart_);
      if (criteria_->className() != "CriteriaMayer") {
        pair_->update(mpart_, 0, "update");
      } else {
        pair_->updatePeTot(pe);
      }
      if (verbose_ == 1) cout << "accepted " << de_ << std::endl;
      trialAccept_();
    } else {
      space()->restore(mpart_);
      if (space()->cellType() > 0) {
        space()->updateCellofiMol(space()->mol()[mpart_.front()]);
      }
      if (verbose_ == 1) cout << "rejected " << de_ << std::endl;
      trialReject_();
    }

  // if preFac is zero, still call criteria for WLTMMC to store old state
  } else {
    if (verbose_ == 1) {
      cout << "move rejected because 0 prefactor " << de_ << endl;
    }
    reject_ = 1;
    lnpMet_ = std::numeric_limits<double>::min();
    de_ = 0;
    criteria_->accept(lnpMet_, pair_->peTot() + de_,
                      trialType_.c_str(), reject_);
    trialReject_();
  }
}

void Trial::updateMaxMoveParam_(const double percent,
  const double upperLimit,
  const double lowerLimit,
  const double targAcceptPer) {
  if (acceptPer() < targAcceptPer) {
    if (percent > 0) {
      if (maxMoveParam > lowerLimit/(1 - percent)) {
        maxMoveParam *= 1 - percent;
      }
    } else {
      if (maxMoveParam < upperLimit/(1 - percent)) {
        maxMoveParam *= 1 - percent;
      }
    }
  }
  if (acceptPer() > targAcceptPer) {
    if (percent > 0) {
      if (maxMoveParam < upperLimit/(1 - percent)) {
        maxMoveParam *= 1 + percent;
      }
    } else {
      if (maxMoveParam > lowerLimit/(1 - percent)) {
        maxMoveParam *= 1 + percent;
      }
    }
  }
  zeroStat();
}


string Trial::printStat(const bool header) {
  stringstream stat;
  if (header) {
    stat << className_ << " ";
  } else {
    stat << acceptPer() << " ";
  }
  if (maxMoveFlag == 1) {
    if (header) {
      stat << "maxMove ";
    } else {
      stat << maxMoveParam << " ";
    }
  }
  return stat.str();
}

void Trial::confine(const double upper, const double lower,
  const int dimension) {
  ASSERT(upper > lower, "upper(" << upper << ") must be greater than lower("
    << lower << ")");
  ASSERT(dimension < space()->dimen(), "given dimension(" << dimension << ") must be "
    << "lower than spatial dimension(" << space()->dimen() << ")");
  confineFlag_ = 1;
  confineUpper_ = upper;
  confineLower_ = lower;
  confineDim_ = dimension;
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_



