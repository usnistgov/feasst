#include "./trial.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

Trial::Trial() {
  defaultConstruction();
}
Trial::Trial(Space *space,
             Pair *pair,
             Criteria *criteria)
  : space_(space),
    pair_(pair),
    criteria_(criteria) {
  defaultConstruction();
}
Trial::Trial(Space *space,
             Pair *pair,
             Criteria *criteria,
             const char* fileName)
  : space_(space),
    pair_(pair),
    criteria_(criteria) {
  ASSERT(fileExists(fileName),
         "restart file(" << fileName << ") doesn't exist");
  defaultConstruction();

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
}

/**
 * defaults in constructor
 */
void Trial::defaultConstruction() {
  className_.assign("Trial");
  zeroStat();
  maxMoveParam = -1;
  maxMoveFlag = 0;
  verbose_ = 0;
  nf_ = 0;
  rAbove_ = -1;
  rBelow_ = -1;
  avbOn_ = false;
}

/**
 * reset object pointers
 */
void Trial::reconstruct(Space* space, Pair *pair, Criteria *criteria) {
  space_ = space;
  pair_ = pair;
  criteria_ = criteria;
  Base::reconstruct();
}

/**
 * Initialize trial move counters for statistics
 */
void Trial::zeroStat() {
  de_ = 0;
  deTot_ = 0;
  accepted_ = 0;
  attempted_ = 0;
}

/**
 * call when attempting a trial
 */
void Trial::attempt() {
  ++attempted_;
  criteria_->store(space_, pair_);
  preFac_ = 0;
  lnpMet_ = 0;
  reject_ = 0;
  de_ = 0.;
  def_ = 0.;
  attempt1();
}

/**
 * call when accepting a trial
 */
void Trial::trialAccept() {
  ++accepted_;
  deTot_ += de_;
  WARN(verbose_ == 1, "accepted " << de_);
}

/**
 * call when rejecting a trial
 */
void Trial::trialReject() {
  WARN(verbose_ == 1, "rejected " << de_);
  de_ = 0;
}

/**
 * set the number of first bead insertion attempts
 */
void Trial::numFirstBeads(const int nf    //!< number of first bead attempts
  ) {
  ASSERT((nf == 0) || (nf > 1),
         "nf(" << nf << ") should be 0 or greater than 1");
  nf_ = nf;
  en_.resize(nf_);
  w_.resize(nf_);
  cpdf_.resize(nf_);
}

/**
 * initialize aggregation volume bias
 */
void Trial::initAVB(const double rAbove,   //!< upper limit of bond
                    const double rBelow    //!< lower limit of bond
  ) {
  ASSERT(rAbove > rBelow, "Null aggregation volume when rAbove = "
         << rAbove << " and rBelow = " << rBelow);
  ASSERT(rAbove <= 0.5*space_->minl(), "aggregation volume upper radius("
    << rAbove << ") cannot extend beyond periodic boundary conditions when"
    << "minimum box length is " << space_->minl());
  ASSERT(rAbove <= pair_->rCut(), "aggregation volume upper radius("
    << rAbove << ") cannot extend beyond pair cut-off " << pair_->rCut());
  rAbove_ = rAbove;
  rBelow_ = rBelow;
  vIn_ = volShell(rAbove_, rBelow_, space_->dimen());
  avbOn_ = true;
}

/**
 * compute rosenbluth weights of multiple first bead attempts
 *  and select
 */
double Trial::multiFirstBead(const int flag) {
  ASSERT((flag == 0) || (flag == 1) || (flag == 2) || (flag == 3),
    "Unrecognized flag for multiple first beads");

  // record position
  space_->xStoreMulti(mpart_, -1);
  pair_->cheapEnergy(1);
  for (int i = 0; i < nf_; ++i) {
    if (i != 0) {
      if (avbOn_) {
        space_->avb(mpart_, tmpart_, rAbove_, rBelow_, region_.c_str());
      } else {
        // -1 flag for displacement by half box length in any direction
        space_->randDisp(mpart_, -1);
        // -1 flag for completely random displacement
        space_->randRotate(mpart_, -1);
      }
      space_->wrap(mpart_);
      // -2 flag to store multiple positions
      space_->xStoreMulti(mpart_, -2);
    }
    en_[i] = pair_->multiPartEner(mpart_, flag);
    w_[i] = exp(-criteria_->beta()*en_[i]);
    if (i != 0) w_[i] += w_[i-1];
  }
  const double wr = w_.back();
  // cout << "wr " << wr << endl;
  pair_->cheapEnergy(0);

  // select trial position from rosenblueth factor weights
  if (wr != 0) {
    // if new configuration, select f and multiply by weight
    if ( (flag == 1) || (flag == 3) ) {
      for (int i = 0; i < nf_; ++i) cpdf_[i] = w_[i]/wr;
      int f = ranFromCPDF(cpdf_);
      def_ = en_[f];
      space_->xStoreMulti(mpart_, f);
      return wr/static_cast<double>(nf_);

    // if old configuration, select first trial and divide by weight
    } else {
      def_ = en_[0];
      space_->xStoreMulti(mpart_, 0);
      // cout << "def_ " << def_ << " nf/wr " << nf_/wr << endl;
      return static_cast<double>(nf_)/wr;
    }
  } else {
    def_ = 0;
    return 0;
  }
}

/**
 * replace or restore criteria pointer
 */
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

/**
 * write restart file
 */
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
  if (maxMoveParam != -1) {
    file << "# maxMoveParam " << maxMoveParam << endl;
  }

  // write random number generator state
  writeRngRestart(fileName);
}

/**
 * return acceptance percentage
 */
double Trial::acceptPer() const {
  if (attempted_ != 0) {
    return static_cast<double>(accepted_)/attempted_;
  } else {
    return 0;
  }
}

/**
 * call when recording old configuration
 */
void Trial::trialMoveRecord() {
  peOld_ = pair_->multiPartEner(mpart_, 0);
  pair_->update(mpart_, 0, "store");
  space_->xStore(mpart_);
}

/**
 * call when recording old configuration before collective move
 */
void Trial::trialMoveRecordAll(const int flag) {
  peOld_ = pair_->allPartEnerForce(flag);
  pair_->update(space_->listAtoms(), 0, "store");
  space_->xStoreAll();
}

/**
 * call to decide whether to accept or reject trial move
 */
void Trial::trialMoveDecide(const double def,   //!< energy of first bead
  const double preFac   //!< acceptance criteria prefactor
  ) {
  if (preFac != 0) {
    // compute energy contribution of selected molecule in new configuration
    if (space_->cellType() > 0) {
      space_->updateCellofiMol(space_->mol()[mpart_.front()]);
    }
    const double pe = pair_->multiPartEner(mpart_, 1);
    de_ = pe - peOld_;
    lnpMet_ = log(preFac) - criteria_->beta()*(de_ - def);
    reject_ = 0;
    // if accepted, update
    if (criteria_->accept(lnpMet_, pair_->peTot() + de_, trialType_.c_str(),
                          reject_) == 1) {
      space_->wrap(mpart_);
      pair_->update(mpart_, 0, "update");
      if (verbose_ == 1) cout << "accepted " << de_ << std::endl;
      trialAccept();
    } else {
      space_->restore(mpart_);
      if (space_->cellType() > 0) {
        space_->updateCellofiMol(space_->mol()[mpart_.front()]);
      }
      if (verbose_ == 1) cout << "rejected " << de_ << std::endl;
      trialReject();
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
    trialReject();
  }
}

/**
 * update the maxMoveParam based on percentage change and limits
 */
void Trial::updateMaxMoveParam(const double percent,
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


/*
 * return string for status of trial
 */
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

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_



