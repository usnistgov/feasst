#include "./trial_transform.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

TrialTransform::TrialTransform(
  const char* transType)    //!< type of transformation
  : Trial(),
    transType_(transType) {
  defaultConstruction();
}
TrialTransform::TrialTransform(Space *space,
  Pair *pair,
  Criteria *criteria,
  const char* transType)    //!< type of transformation
  : Trial(space, pair, criteria),
    transType_(transType) {
  defaultConstruction();
}
TrialTransform::TrialTransform(const char* fileName,
  Space *space,
  Pair *pair,
  Criteria *criteria)
  : Trial(space, pair, criteria, fileName) {
  transType_ = fstos("transType", fileName);
  defaultConstruction();
  targAcceptPer = fstod("targAcceptPer", fileName);

  // although maxMoveParam was already read in the base class
  // read it again because it was over-written by defaultConstruction
  maxMoveParam = fstod("maxMoveParam", fileName);
}

/**
 * write restart file
 */
void TrialTransform::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  file << "# transType " << transType_ << endl;
  file << "# targAcceptPer " << targAcceptPer << endl;
}

/**
 * default construction
 */
void TrialTransform::defaultConstruction() {
  className_.assign("TrialTransform");
  trialType_.assign("move");
  verbose_ = 0;
  maxMoveFlag = 1;
  if ( (transType_.compare("translate") == 0) ||
       (transType_.compare("rotate") == 0) ) {
    maxMoveParam = 0.1;
    targAcceptPer = 0.25;
  } else if (transType_.compare("smctrans") == 0) {
    maxMoveParam = 2e-4;
    targAcceptPer = 0.25;
  } else if ( (transType_.compare("xytilt") == 0) ||
              (transType_.compare("xztilt") == 0) ||
              (transType_.compare("yztilt") == 0) ) {
    maxMoveParam = 0.1;
    targAcceptPer = 0.25;
  } else if ( (transType_.compare("vol") == 0) ||
              (transType_.compare("lxmod") == 0) ||
              (transType_.compare("lymod") == 0) ||
              (transType_.compare("lzmod") == 0) ) {
    maxMoveParam = 0.1;
    targAcceptPer = 0.25;
  } else {
    ASSERT(0, "transformation type (" << transType_ << ") not recognized");
  }
}

/**
 * Attempt trial
 */
void TrialTransform::attempt1() {
  if (verbose_ == 1) {
    cout << std::setprecision(std::numeric_limits<double>::digits10+2)
         << "attempting " << transType_ << " " << pair_->peTot() << endl;
  }
  if (space_->nMol() > 0) {
    if (transType_.compare("smctrans") == 0) {
      trialMoveRecordAll(1);

      // rigidily translate molecules according to the force on their center
      // of mass
      const double A = maxMoveParam*maxMoveParam/2.;
      const double betaA = criteria_->beta()*A;
      vector<vector<double> > dr(space_->nMol(),
                                 vector<double>(space_->dimen()));
      double wold = 0;    //!< prefactor for old configuration
      vector<double> R;   //!< multivariate guassian random numbers
      const vector<vector<double> > &fCOM = pair_->fCOM();
      double dx;
      const vector<double> l = space_->l();
      bool tooBig = false;
      for (int iMol = 0; iMol < space_->nMol(); ++iMol) {
        // generate displacement vector from fCOM
        R = stdRanNum(maxMoveParam, space_->dimen());
        for (int dim = 0; dim < space_->dimen(); ++dim) {
          dr[iMol][dim] = dx = betaA*fCOM[iMol][dim] + R[dim];
          if (fabs(dx) > 0.5 * l[dim]) tooBig = true;
        }

        // translate molecule
        space_->transMol(iMol, dr[iMol]);

        // compute bias prefactor from old forces
        wold -= vecDotProd(R, R);
      }
      wold = exp(wold/(4.*A));
      if ( (tooBig == true) || (wold == 0) ) {
        reject_ = 1;
        de_ = 0;
        lnpMet_ = std::numeric_limits<double>::min();
      } else {
        // compute forces and energy of new configuration
        if (space_->cellType() > 0) space_->updateCellofallMol();
        de_ = pair_->allPartEnerForce(1) - peOld_;

        // compute bias prefactor from new forces
        double wnew = 0;
        for (int iMol = 0; iMol < space_->nMol(); ++iMol) {
          for (int dim = 0; dim < space_->dimen(); ++dim) {
            R[dim] = -dr[iMol][dim] - betaA*fCOM[iMol][dim];
          }
          wnew -= vecDotProd(R, R);
        }
        wnew = exp(wnew/(4.*A));
        lnpMet_ = log(wnew/wold) - criteria_->beta()*de_;
        reject_ = 0;
      }

      // accept or reject with bias prefactor
      if (criteria_->accept(lnpMet_, pair_->peTot() + de_, trialType_.c_str(),
          reject_) == 1) {
        for (int iMol = 0; iMol < space_->nMol(); ++iMol) {
          space_->wrap(space_->imol2mpart(iMol));
        }
        if (pair_->neighOn()) pair_->buildNeighList();
        pair_->update(space_->listAtoms(), 0, "update");
        // cout << "accepted " << transType_ << " " << de_ << endl;
        trialAccept();
      } else {
        space_->restoreAll();
        if (space_->cellType() > 0) space_->updateCellofallMol();
        // cout << "rejected " << transType_ << " " << de_ << endl;
        trialReject();
      }

    // floppy box
    } else if ( (transType_.compare("xytilt") == 0) ||
                (transType_.compare("xztilt") == 0) ||
                (transType_.compare("yztilt") == 0) ) {
      trialMoveRecordAll(0);

      // randomly attempt to increase or decrease tilt by maxMoveParam
      //  this must be accompanied by a transformation of the particles
      double dxyt = maxMoveParam*(2*uniformRanNum()-1);
      double tiltOld = 0.;
      if (transType_.compare("xytilt") == 0) {
        tiltOld = space_->xyTilt();
        space_->modXYTilt(dxyt);
      } else if (transType_.compare("xztilt") == 0) {
        tiltOld = space_->xzTilt();
        space_->modXZTilt(dxyt);
      } else if (transType_.compare("yztilt") == 0) {
        tiltOld = space_->yzTilt();
        space_->modYZTilt(dxyt);
      }

      // compute energy of new configuration
      de_ = pair_->allPartEnerForce(1) - peOld_;
      lnpMet_ = -criteria_->beta()*de_;

      // accept or reject with bias prefactor
      if (criteria_->accept(lnpMet_, pair_->peTot() + de_, trialType_.c_str(),
          reject_) == 1) {
        if (pair_->neighOn()) pair_->buildNeighList();
        pair_->update(space_->listAtoms(), 0, "update");
        space_->wrapMol();
        trialAccept();
      } else {
        space_->restoreAll();
        if (transType_.compare("xytilt") == 0) {
          space_->setXYTilt(tiltOld);
        } else if (transType_.compare("xztilt") == 0) {
          space_->setXZTilt(tiltOld);
        } else if (transType_.compare("yztilt") == 0) {
          space_->setYZTilt(tiltOld);
        }
        ASSERT(space_->cellType() <= 0,
          "xytilt trial move not implemented correctly with cell list");
        if (space_->cellType() > 0) space_->updateCellofallMol();
        // cout << "rejected " << transType_ << " " << de_ << endl;
        trialReject();
      }

    // box size change
    } else if ( (transType_.compare("vol") == 0) ||
                (transType_.compare("lxmod") == 0) ||
                (transType_.compare("lymod") == 0) ||
                (transType_.compare("lzmod") == 0) ) {
      trialMoveRecordAll(0);

      // randomly attempt to increase or decrease (lx,ly,lz) by maxMoveParam
      //  this must be accompanied by a transformation of the particles
      const double dlnv = maxMoveParam*(2*uniformRanNum()-1),
       vOld = space_->vol(),
       fac = exp(log(vOld) + dlnv)/vOld;
      scaleAttempt_(fac);
      space_->wrapMol();

      // if box attempts to go beyond certain bounds, it may be modified
      const double facActual = space_->vol()/vOld;

      // compute energy of new configuration
      de_ = pair_->allPartEnerForce(1) - peOld_;
      lnpMet_ = -criteria_->beta()*(de_ + criteria_->pressure()*(vOld*
        (facActual-1.)) - (space_->nMol()+1)*log(facActual)/criteria_->beta());
      // cout << "de " << de_ << "pmet " << lnpMet_ << endl;

      // accept or reject with bias prefactor
      if (criteria_->accept(lnpMet_, pair_->peTot() + de_, trialType_.c_str(),
                            reject_) == 1) {
        if (pair_->neighOn()) pair_->buildNeighList();
        pair_->update(space_->listAtoms(), 0, "update");
        trialAccept();
      } else {
        scaleAttempt_(1./facActual);
        space_->restoreAll();
        ASSERT(space_->cellType() <= 0,
          "xytilt trial move not implemented correctly with cell list");
        if (space_->cellType() > 0) space_->updateCellofallMol();
        // cout << "rejected " << transType_ << " " << de_ << endl;
        trialReject();
      }

      // record statistics
      if (reject_ != 1) {
        if (transType_.compare("vol") == 0) {
          paramAccumulator.accumulate(space_->vol());
        } else if (transType_.compare("lxmod") == 0) {
          paramAccumulator.accumulate(space_->l(0));
        } else if (transType_.compare("lymod") == 0) {
          paramAccumulator.accumulate(space_->l(1));
        } else if (transType_.compare("lzmod") == 0) {
          paramAccumulator.accumulate(space_->l(2));
        }
      }

    // rigid translation or rotation
    } else {
      mpart_ = space_->randMol();    // select a random molecule
      trialMoveRecord();
      if (transType_.compare("translate") == 0) {
        space_->randDisp(mpart_, maxMoveParam);
      } else if (transType_.compare("rotate") == 0) {
        space_->randRotate(mpart_, maxMoveParam);
      }
      trialMoveDecide(0, 1);
    }
  } else {
    // ensured rejection, however, criteria can update
    trialMoveDecide(0, 0);
  }
}

/**
 * tune max parameters for translation and rotation moves
 */
void TrialTransform::tuneParameters() {
  // determine limits and percentage changes
  double percent = 0.05;   // percentage change
  double upperLimit = 0;         // do not change above this limit
  double lowerLimit = 1e-5;
  if ( (transType_.compare("translate") == 0) ||
       (transType_.compare("smctrans") == 0) ||
       (transType_.compare("xytilt") == 0) ||
       (transType_.compare("xztilt") == 0) ||
       (transType_.compare("yztilt") == 0) ||
       (transType_.compare("vol") == 0) ||
       (transType_.compare("lxmod") == 0) ||
       (transType_.compare("lymod") == 0) ||
       (transType_.compare("lzmod") == 0)
     ) {
    upperLimit = space_->minl()/4.;
  } else if ( (transType_.compare("rotate") == 0) ) {
    upperLimit = 1e1;
  } else {
    ASSERT(0, "unrecognized transType(" << transType_ << ")");
  }

  if (upperLimit != 0) {
    updateMaxMoveParam(percent, upperLimit, lowerLimit, targAcceptPer);
  }
}

/*
 * return string for status of trial
 */
string TrialTransform::printStat(const bool header) {
  stringstream stat;
  if (header) {
    stat << transType_ << " ";
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
  if (transType_.compare("vol") == 0) {
    if (header) {
      stat << "volume ";
    } else {
      stat << space_->vol() << " ";
    }
  }
  if (transType_.compare("lxmod") == 0) {
    if (header) {
      stat << "lx ";
    } else {
      stat << space_->l(0) << " ";
    }
  }
  if (transType_.compare("lymod") == 0) {
    if (header) {
      stat << "ly ";
    } else {
      stat << space_->l(1) << " ";
    }
  }
  if (transType_.compare("lzmod") == 0) {
    if (header) {
      stat << "lz ";
    } else {
      stat << space_->l(2) << " ";
    }
  }
  return stat.str();
}

/*
 * attempt to scale the domain
 */
void TrialTransform::scaleAttempt_(const double factor) {
  // determine if lx, ly or lz (or volume if dim == -1)
  if (transType_.compare("lxmod") == 0) {
    space_->scaleDomain(factor, 0);
  } else if (transType_.compare("lymod") == 0) {
    space_->scaleDomain(factor, 1);
  } else if (transType_.compare("lzmod") == 0) {
    space_->scaleDomain(factor, 2);
  } else if (transType_.compare("vol") == 0) {
    space_->scaleDomain(factor);
  } else {
    ASSERT(0, "unrecognized transType_(" << transType_ << ") for scale.");
  }
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_



