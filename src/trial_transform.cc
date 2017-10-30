/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./trial_transform.h"
#include "mc.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

TrialTransform::TrialTransform(
  const char* transType)
  : Trial(),
    transType_(transType) {
  defaultConstruction_();
}

TrialTransform::TrialTransform(
  Pair *pair,
  Criteria *criteria,
  const char* transType)
  : Trial(pair, criteria),
    transType_(transType) {
  defaultConstruction_();
}

TrialTransform::TrialTransform(const char* fileName,
  Pair *pair,
  Criteria *criteria)
  : Trial(pair, criteria, fileName) {
  transType_ = fstos("transType", fileName);
  defaultConstruction_();
  targAcceptPer = fstod("targAcceptPer", fileName);

  // although maxMoveParam was already read in the base class
  // read it again because it was over-written by defaultConstruction
  maxMoveParam = fstod("maxMoveParam", fileName);

  string strtmp = fstos("molType", fileName);
  if (!strtmp.empty()) {
    molType_ = strtmp;
  }
}

void TrialTransform::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  file << "# transType " << transType_ << endl;
  file << "# targAcceptPer " << targAcceptPer << endl;
  if (!molType_.empty()) {
    file << "# molType " << molType_ << endl;
  }
}

void TrialTransform::defaultConstruction_() {
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

void TrialTransform::attempt1_() {
  if (verbose_ == 1) {
    cout << std::setprecision(std::numeric_limits<double>::digits10+2)
         << "attempting " << transType_ << " " << pair_->peTot() << endl;
  }
  if (space()->nMol() > 0) {
    if (transType_.compare("smctrans") == 0) {
      trialMoveRecordAll_(1);

      // rigidily translate molecules according to the force on their center
      // of mass
      const double A = maxMoveParam*maxMoveParam/2.;
      const double betaA = criteria_->beta()*A;
      vector<vector<double> > dr(space()->nMol(),
                                 vector<double>(space()->dimen()));
      double wold = 0;    //!< prefactor for old configuration
      vector<double> R;   //!< multivariate guassian random numbers
      const vector<vector<double> > &fCOM = pair_->fCOM();
      double dx;
      const vector<double> l = space()->l();
      bool tooBig = false;
      for (int iMol = 0; iMol < space()->nMol(); ++iMol) {
        // generate displacement vector from fCOM
        R = stdRanNum(maxMoveParam, space()->dimen());
        for (int dim = 0; dim < space()->dimen(); ++dim) {
          dr[iMol][dim] = dx = betaA*fCOM[iMol][dim] + R[dim];
          if (fabs(dx) > 0.5 * l[dim]) tooBig = true;
        }

        // translate molecule
        space()->transMol(iMol, dr[iMol]);

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
        if (space()->cellType() > 0) space()->updateCellofallMol();
        de_ = pair_->allPartEnerForce(1) - peOld_;

        // compute bias prefactor from new forces
        double wnew = 0;
        for (int iMol = 0; iMol < space()->nMol(); ++iMol) {
          for (int dim = 0; dim < space()->dimen(); ++dim) {
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
        for (int iMol = 0; iMol < space()->nMol(); ++iMol) {
          space()->wrap(space()->imol2mpart(iMol));
        }
        if (pair_->neighOn()) pair_->buildNeighList();
        pair_->update(space()->listAtoms(), 0, "update");
        // cout << "accepted " << transType_ << " " << de_ << endl;
        trialAccept_();
      } else {
        space()->restoreAll();
        if (space()->cellType() > 0) space()->updateCellofallMol();
        // cout << "rejected " << transType_ << " " << de_ << endl;
        trialReject_();
      }

    // floppy box
    } else if ( (transType_.compare("xytilt") == 0) ||
                (transType_.compare("xztilt") == 0) ||
                (transType_.compare("yztilt") == 0) ) {
      trialMoveRecordAll_(0);

      // randomly attempt to increase or decrease tilt by maxMoveParam
      //  this must be accompanied by a transformation of the particles
      double dxyt = maxMoveParam*(2*uniformRanNum()-1);
      double tiltOld = 0.;
      if (transType_.compare("xytilt") == 0) {
        tiltOld = space()->xyTilt();
        space()->modXYTilt(dxyt);
      } else if (transType_.compare("xztilt") == 0) {
        tiltOld = space()->xzTilt();
        space()->modXZTilt(dxyt);
      } else if (transType_.compare("yztilt") == 0) {
        tiltOld = space()->yzTilt();
        space()->modYZTilt(dxyt);
      }

      // compute energy of new configuration
      de_ = pair_->allPartEnerForce(1) - peOld_;
      lnpMet_ = -criteria_->beta()*de_;

      // accept or reject with bias prefactor
      if (criteria_->accept(lnpMet_, pair_->peTot() + de_, trialType_.c_str(),
          reject_) == 1) {
        if (pair_->neighOn()) pair_->buildNeighList();
        pair_->update(space()->listAtoms(), 0, "update");
        space()->wrapMol();
        trialAccept_();
      } else {
        space()->restoreAll();
        if (transType_.compare("xytilt") == 0) {
          space()->setXYTilt(tiltOld);
        } else if (transType_.compare("xztilt") == 0) {
          space()->setXZTilt(tiltOld);
        } else if (transType_.compare("yztilt") == 0) {
          space()->setYZTilt(tiltOld);
        }
        ASSERT(space()->cellType() <= 0,
          "xytilt trial move not implemented correctly with cell list");
        if (space()->cellType() > 0) space()->updateCellofallMol();
        // cout << "rejected " << transType_ << " " << de_ << endl;
        trialReject_();
      }

    // box size change
    } else if ( (transType_.compare("vol") == 0) ||
                (transType_.compare("lxmod") == 0) ||
                (transType_.compare("lymod") == 0) ||
                (transType_.compare("lzmod") == 0) ) {
      trialMoveRecordAll_(0);

      // randomly attempt to increase or decrease (lx,ly,lz) by maxMoveParam
      //  this must be accompanied by a transformation of the particles
      const double dlnv = maxMoveParam*(2*uniformRanNum()-1),
       vOld = space()->vol(),
       fac = exp(log(vOld) + dlnv)/vOld;
//      cout << "fac " << fac << " vOld " << vOld << " pres " << criteria_->pressure() << " " << space()->xyTilt() << " " << space()->xzTilt() << " " << space()->yzTilt() << " " << space()->l(0) << " " << space()->l(1) << " " << space()->l(2) << endl;
      scaleAttempt_(fac);
      space()->wrapMol();

      // if box attempts to go beyond certain bounds, it may be modified
      const double facActual = space()->vol()/vOld;

      // compute energy of new configuration
      de_ = pair_->allPartEnerForce(1) - peOld_;
      lnpMet_ = -criteria_->beta()*(de_ + criteria_->pressure()*(vOld*
        (facActual-1.)) - (space()->nMol()+1)*log(facActual)/criteria_->beta());
      // cout << "de " << de_ << "pmet " << lnpMet_ << endl;

      // accept or reject with bias prefactor
      if (criteria_->accept(lnpMet_, pair_->peTot() + de_, trialType_.c_str(),
                            reject_) == 1) {
        if (pair_->neighOn()) pair_->buildNeighList();
        pair_->update(space()->listAtoms(), 0, "update");
        trialAccept_();
      } else {
        scaleAttempt_(1./facActual);
        space()->restoreAll();
        ASSERT(space()->cellType() <= 0,
          "xytilt trial move not implemented correctly with cell list");
        if (space()->cellType() > 0) space()->updateCellofallMol();
        // cout << "rejected " << transType_ << " " << de_ << endl;
        trialReject_();
      }

      // record statistics
      if (reject_ != 1) {
        if (transType_.compare("vol") == 0) {
          paramAccumulator.accumulate(space()->vol());
        } else if (transType_.compare("lxmod") == 0) {
          paramAccumulator.accumulate(space()->l(0));
        } else if (transType_.compare("lymod") == 0) {
          paramAccumulator.accumulate(space()->l(1));
        } else if (transType_.compare("lzmod") == 0) {
          paramAccumulator.accumulate(space()->l(2));
        }
      }

    // rigid translation or rotation
    } else {
      if (molType_.empty()) {
        // if no molType given, move any particle
        mpart_ = space()->randMol();    // select a random molecule
      } else {
        // otherwise, move only molType particles
        const int iMolIndex = space()->findAddMolListIndex(molType_);
        const int nMolOfType = space()->nMolType()[iMolIndex];
        if (nMolOfType > 0) {
          const int iMol = space()->randMolofType(iMolIndex);
          if (iMol == -1) {
            reject_ = 1;
          } else {
            mpart_ = space()->imol2mpart(iMol);
          }
        } else {
          reject_ = 1;
        }
      }
      if (reject_ == 1) {
        // ensured rejection, however, criteria can update
        trialMoveDecide_(0, 0);
      } else {
        trialMoveRecord_();
        if (transType_.compare("translate") == 0) {
          space()->randDisp(mpart_, maxMoveParam);
        } else if (transType_.compare("rotate") == 0) {
          space()->randRotate(mpart_, maxMoveParam);
        }
        trialMoveDecide_(0, 1);
      }
    }
  } else {
    // ensured rejection, however, criteria can update
    trialMoveDecide_(0, 0);
  }
}

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
    upperLimit = space()->minl()/4.;
  } else if ( (transType_.compare("rotate") == 0) ) {
    upperLimit = 1e1;
  } else {
    ASSERT(0, "unrecognized transType(" << transType_ << ")");
  }

  if (upperLimit != 0) {
    updateMaxMoveParam_(percent, upperLimit, lowerLimit, targAcceptPer);
  }
}

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
      stat << space()->vol() << " ";
    }
  }
  if (transType_.compare("lxmod") == 0) {
    if (header) {
      stat << "lx ";
    } else {
      stat << space()->l(0) << " ";
    }
  }
  if (transType_.compare("lymod") == 0) {
    if (header) {
      stat << "ly ";
    } else {
      stat << space()->l(1) << " ";
    }
  }
  if (transType_.compare("lzmod") == 0) {
    if (header) {
      stat << "lz ";
    } else {
      stat << space()->l(2) << " ";
    }
  }
  return stat.str();
}

void TrialTransform::scaleAttempt_(const double factor) {
  ASSERT( (factor > DTOL) && (factor < NUM_INF),
    "cannot scale domain by a factor: " << factor);
  // determine if lx, ly or lz (or volume if dim == -1)
  if (transType_.compare("lxmod") == 0) {
    space()->scaleDomain(factor, 0);
  } else if (transType_.compare("lymod") == 0) {
    space()->scaleDomain(factor, 1);
  } else if (transType_.compare("lzmod") == 0) {
    space()->scaleDomain(factor, 2);
  } else if (transType_.compare("vol") == 0) {
    space()->scaleDomain(factor);
  } else {
    ASSERT(0, "unrecognized transType_(" << transType_ << ") for scale.");
  }
}

shared_ptr<TrialTransform> makeTrialTransform(Pair *pair,
  Criteria *criteria, const char* transType) {
  return make_shared<TrialTransform>(pair, criteria, transType);
}

shared_ptr<TrialTransform> makeTrialTransform(const char* transType) {
  return make_shared<TrialTransform>(transType);
}

void transformTrial(MC *mc, const char* type, double maxMoveParam) {
  shared_ptr<TrialTransform> trial = make_shared<TrialTransform>(type);
  if (maxMoveParam != -1) trial->maxMoveParam = maxMoveParam;
  mc->initTrial(trial);
}
void transformTrial(shared_ptr<MC> mc, const char* type,
  double maxMoveParam) {
  transformTrial(mc.get(), type, maxMoveParam);
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_



