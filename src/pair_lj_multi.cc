/**
 * \file
 *
 * \brief lennard-jones pairwise interactions
 *
 * Slow implementation of pair_lj class
 */

#include "./pair_lj_multi.h"

namespace feasst {

/**
 * Constructor for pair_lj class requires the following
 */
PairLJMulti::PairLJMulti(Space* space,
         const double rCut  //!< interaction cut-off distance
  )
    : PairLJ(space, rCut) {
  defaultConstruction();
}
PairLJMulti::PairLJMulti(Space* space,
         const char* fileName
  )
    : PairLJ(space, fileName) {
  defaultConstruction();

  // set exponential type
  string str = fstos("expType", fileName);
  if (!str.empty()) {
    initExpType(stoi(str));
  } else {
    initExpType(0);
  }

  // set yukawa
  str = fstos("yukawaFlag", fileName);
  if (!str.empty()) {
    yukawa_ = stoi(str);
    yukawaA_ = fstoull("yukawaA", fileName);
    yukawaK_ = fstoull("yukawaK", fileName);
  } else {
    yukawa_ = 0;
    yukawaA_ = 0;
    yukawaK_ = 0;
  }

  for (unsigned int i = 0; i < epsij_.size(); ++i) {
    for (unsigned int j = i; j < epsij_.size(); ++j) {
      stringstream ss;
      ss << "peLinearShift" << i << j;
      string strtmp = fstos(ss.str().c_str(), fileName);
      if (!strtmp.empty()) {
        linearShiftijset(i, j, 1);
      } else {
        ss.str("");
        ss << "peShift" << i << j;
        strtmp = fstos(ss.str().c_str(), fileName);
        if (!strtmp.empty()) {
          cutShiftijset(i, j, 1);
        }
      }
    }
  }

  str = fstos("nGaussianParams", fileName);
  gaussian_ = 0;
  if (!str.empty()) {
    gaussian_ = 1;
    const int n = stoi(str);
    for (int i = 0; i < n; ++i) {
      vector<double> params;
      stringstream ss;
      ss << "gausParam" << i << "height";
      params.push_back(fstod(ss.str().c_str(), fileName));
      ss.str("");
      ss << "gausParam" << i << "position";
      params.push_back(fstod(ss.str().c_str(), fileName));
      ss.str("");
      ss << "gausParam" << i << "spread";
      params.push_back(fstod(ss.str().c_str(), fileName));
      gausParam_.push_back(params);
    }
  }

  str = fstos("lambdaFlag", fileName);
  if (!str.empty()) {
    lambdaFlag_ = stoi(str);
    lambda_.resize(space_->nParticleTypes(), vector<double>(
      space_->nParticleTypes()));
    for (int i = 0; i < space_->nParticleTypes(); ++i) {
      for (int j = 0; j < space_->nParticleTypes(); ++j) {
        stringstream ss;
        ss << "lambdaParami" << i << "j" << j;
        str = fstos(ss.str().c_str(), fileName);
        lambda_[i][j] = stod(str);
      }
    }
  } else {
    lambdaFlag_ = 0;
  }

  str = fstos("cutShiftFlag", fileName);
  if (!str.empty()) cutShift(stoi(str));
  if (cutShiftFlag_ == 0) cutShift(0);
}

/**
 * defeaults in constructor
 */
void PairLJMulti::defaultConstruction() {
  className_.assign("PairLJMulti");
  initExpType(0);
  yukawa_ = 0;
  yukawaA_ = 0;
  yukawaK_ = 0;
  lambdaFlag_ = 0;
  gaussian_ = 0;
}

/**
 * Lennard-Jones pair-wise force calculation
 */
int PairLJMulti::initEnergy() {
  const int verbose = 0;
  if (verbose) cout << "in pair_lj_multi.cc" << endl;

  // shorthand for read-only space variables
  const int dimen = space_->dimen();
  const int natom = space_->natom();
  const vector<double> l = space_->l();
  const vector<double> x = space_->x();
  const vector<int> mol = space_->mol();
  const vector<int> type = space_->type();

  // zero accumulators: potential energy, force, and virial
  std::fill(pe_.begin(), pe_.end(), 0.);
  myFill(0., f_);
  myFill(0., vr_);
  peLJ_ = 0;
  fCOM_.clear();
  fCOM_.resize(space_->nMol(), vector<double>(dimen_, 0.));

  double r2inv, r6inv, pePart;

  // loop through nearest neighbor atom pairs
  for (int ipart = 0; ipart < natom - 1; ++ipart) {
    if (eps_[type[ipart]] != 0) {
      for (int jpart = ipart + 1; jpart < natom; ++jpart) {
        if ( intraCheck_(ipart, jpart, mol[ipart], mol[jpart]) &&
             (eps_[type[jpart]] != 0) ) {
          // separation vector, xij with periodic boundary conditions
          vector<double> xij(dimen_);
          for (int i = 0; i < dimen_; ++i) {
            xij[i] = x[dimen_*ipart+i] - x[dimen_*jpart+i];
          }
          const vector<double> dx = space_->pbc(xij);
          for (int dim = 0; dim < dimen_; ++dim) {
            xij[dim] += dx[dim];
          }
          double r2 = myVecDotProd(xij, xij);

          if (verbose) cout << "r2 " << r2 << " i " << ipart << " j " << jpart
             << " itype " << type[ipart] << " jtype " << type[jpart] << endl;

          // no interaction beyond cut-off distance
          if (r2 < pow(rCutij_[type[ipart]][type[jpart]], 2.)) {
            // energy
            pePart = 0;
            const double eps = epsij_[type[ipart]][type[jpart]],
                         sig = sigij_[type[ipart]][type[jpart]];
            if (sigrefFlag_ != 1) {
              r2inv = pow(sig, 2)/r2;
            } else {
              const double sigref = sigRefij_[type[ipart]][type[jpart]];
              r2inv = sigref/(sqrt(r2) - sig + sigref);
              r2inv = r2inv*r2inv;

              // inner hard sphere
              if (sqrt(r2) <= sig - sigref) {
                pePart += std::numeric_limits<double>::max()/1e10;
              }
            }
            r6inv = r2inv*r2inv*r2inv;
            if (expType_ == 1) {
              r6inv = r6inv*r6inv;
            } else if (expType_ == 2) {
              r6inv = pow(r2inv, 16.6755*0.5);
            } else if (expType_ == 3) {
              r6inv = pow(r2inv, 25);
            } else if (expType_ == 4) {
              r6inv = pow(r2inv, 64);
            } else if (expType_ == 5) {
              r6inv = pow(r2inv, 12);
            } else if (expType_ == 6) {
              r6inv = pow(r2inv, 9);
            }
            pePart += eps *(4. * (r6inv*(r6inv - 1.))
              + peShiftij_[type[ipart]][type[jpart]]);
            if (linearShiftFlag_) {
              pePart += peLinearShiftij_[type[ipart]][type[jpart]]
                * (sqrt(r2) - rCutij_[type[ipart]][type[jpart]]);
            }
            if (lambdaFlag_ != 0) {
              const double lambda = lambda_[type[ipart]][type[jpart]];
              double rwca;
              if (sigrefFlag_ != 1) {
                rwca = pow(2., 1./alpha_)*sig;
              } else {
                const double sigref = sigRefij_[type[ipart]][type[jpart]];
                rwca = pow(2., 1./alpha_)*sigref + sig - sigref;
              }
              if (r2 < rwca*rwca) {
                pePart += eps*(1. - lambda) + (lambda - 1.)
                  *peShiftij_[type[ipart]][type[jpart]];
                if (linearShiftFlag_) {
                  pePart += (lambda - 1.)*(sqrt(r2)
                    - rCutij_[type[ipart]][type[jpart]])
                    *peLinearShiftij_[type[ipart]][type[jpart]];
                }
              } else {
                pePart *= lambda;
              }
              if (verbose) cout << "pePart " << pePart << endl;
            }

            // yukawa
            if (yukawa_ == 1) {
              const double r = sqrt(r2);
              pePart += eps*yukawaA_ * exp(-yukawaK_*r/sig)/r*sig;
            } else if (yukawa_ == 2) {
              const double r = sqrt(r2);
              pePart += yukawaA_ * exp(-yukawaK_*r/sig)/r*sig;
            }

            // gaussians
            if (gaussian_ == 1) {
              for (unsigned int ig = 0; ig < gausParam_.size(); ++ig) {
                pePart += gausParam_[ig][0]*exp(-pow((sqrt(r2)
                  - gausParam_[ig][1])/gausParam_[ig][2], 2));
              }
            }
            
            peLJ_ += pePart;
            pe_[ipart] += pePart/2.;
            pe_[jpart] += pePart/2.;

            if (verbose) {
              std::streamsize ss = cout.precision();
              cout << std::setprecision(std::numeric_limits<long double>::digits10+2) << "pepart " << pePart << " for ipart " << ipart << " jpart " << jpart << " r6inv " << r6inv << " eps " << epsij_[type[ipart]][type[jpart]] << " peshi " << peShiftij_[type[ipart]][type[jpart]] << " lsf " << linearShiftFlag_ << " sig " << sigij_[type[ipart]][type[jpart]] << " r2 " << r2 << " yuk " << yukawa_ << " yukaA " << yukawaA_ << " yukaK " << yukawaK_ << endl;
              cout << std::setprecision(ss);
            }

            // force
            vector<double> fij(dimen);
            double fPart = 8.*alpha_*epsij_[type[ipart]][type[jpart]]
              * (r6inv*r2inv*(r6inv - 0.5)) / sigij_[type[ipart]][type[jpart]]
              / sigij_[type[ipart]][type[jpart]];
            if (linearShiftFlag_) {
              fPart -= peLinearShiftij_[type[ipart]][type[jpart]]/sqrt(r2);
            }
            if (verbose) cout << "fPart " << fPart << " f " 
              << fPart*sqrt(r2) << endl;
            for (int i = 0; i < dimen; ++i) {
              fij[i] += fPart * xij[i];
              f_[ipart][i] += fij[i];
              f_[jpart][i] -= fij[i];
              fCOM_[mol[ipart]][i] += fij[i];
              fCOM_[mol[jpart]][i] -= fij[i];

              // virial tensor
              for (int j = 0; j < dimen; ++j) {
                vr_[ipart][i][j] += 0.5 * xij[j] * fij[i];
                vr_[jpart][i][j] += 0.5 * xij[j] * fij[i];
              }
            }
          }
        }
      }
    }
  }

  // standard long range corrections
  peLRC_ = 0;
  if (lrcFlag == 1) {
    for (int ipart = 0; ipart < natom; ++ipart) {
      const double enlrc = computeLRC(ipart);
      pe_[ipart] += enlrc;
      peLRC_ += enlrc;
    }
  }

  peTot_ = peLJ_ + peLRC_;
  return 0;
}

/**
 * potential energy and forces of all particles
 *  if flag == 0, dummy calculation
 *  if flag == 1, all config calculation
 */
double PairLJMulti::allPartEnerForce(const int flag) {
  peSRone_ = 0;
  // standard long range corrections
  peLRCone_ = 0;
  if (lrcFlag == 1) {
    for (int ipart = 0; ipart < space_->natom(); ++ipart) {
      const double enlrc = computeLRC(ipart);
      peLRCone_ += enlrc;
    }
  }

  if (flag == 0) {
    peSRone_ = peTot() - peLRCone_;
    return peSRone_ + peLRCone_;
  } else {
    // zero accumulators: potential energy and force
    fCOM_.clear();
    fCOM_.resize(space_->nMol(), vector<double>(dimen_, 0.));

    if ( (space_->cellType() == 0) || (rCutMaxAll_ > space_->dCellMin()) ) {
      if (dimen_ == 3) {
        return allPartEnerForceAtomCutNoCell() + peLRCone_;
      } else if (dimen_ == 2) {
        return allPartEnerForceAtomCutNoCell2D() + peLRCone_;
      }
    } else if (space_->cellType() == 1) {
      return allPartEnerForceAtomCutCell() + peLRCone_;
    } else {
      ASSERT(0, "cell type(" << space_->cellType() << ")");
    }
  }
  return 1e300;
}

/**
 * Lennard-Jones potential energy contribution due to particles
 */
double PairLJMulti::multiPartEner(
  const vector<int> mpart,   //!< particle indx to calculate energy interactions
  const int flag) {     //!< place holder for other pair styles
  if (flag == 0) {}  // remove unused parameter warning

  // zero potential energy contribution of particle ipart
  peSRone_ = 0;
  peLRCone_ = 0;

  // contribution of particles to standard long range corrections
  // (non-additive ~n^2)
  if ( (!cheapEnergy_) && (lrcFlag == 1) ) {
    for (unsigned int ii = 0; ii < mpart.size(); ++ii) {
      const int ipart = mpart[ii];
      const int iType = space_->type()[ipart];

      // search all particle types, and sum lrcs of each
      for (int jType = 0; jType < space_->nParticleTypes(); ++jType) {
        int n = space_->nType()[jType];

// ** NOTE: LRCs include self-term, due to interactions with pbcs
//        // subtract number of jTypes in current molecule
//        const int iMol = space_->mol()[ipart];
//        for (int i = space_->mol2part()[iMol];
//          i < space_->mol2part()[iMol+1]; ++i) {
//          if (space_->type()[i] == jType) --n;
//        }

        peLRCone_ += (2.*static_cast<double>(n)-1)/space_->vol()
          * lrcPreCalc_[iType][jType];
      }
    }
  }

  if (dimen_ == 3) {
    return multiPartEnerAtomCut(mpart) + peLRCone_;
  } else {
    return multiPartEnerAtomCut2D(mpart) + peLRCone_;
  }
}

/**
 * initialize cut and shifted potential
 *  if flag == 0, do not shift
 *  if flag == 1, shift by potential value to zero at rCut
 */
void PairLJMulti::cutShift(const int flag    //!< initialization flag
  ) {
  cutShiftFlag_ = flag;
  peShiftij_.resize(epsij_.size(), vector<double>(epsij_.size()));
  rCutij_.clear();
  rCutij_.resize(epsij_.size(), vector<double>(epsij_.size(), rCut_));
  myFill(rCut_, rCutij_);
  if (flag == 0) {
    myFill(0., peShiftij_);
  } else if (flag == 1) {
    for (unsigned int i = 0; i < epsij_.size(); ++i) {
      for (unsigned int j = 0; j < epsij_.size(); ++j) {
        cutShiftijset(i, j, flag);
      }
    }
    lrcFlag = 0;
  } else {
    ASSERT(0, "Unrecognized cutShift flag(" << flag << ").");
  }
}

/**
 * initialize linear force shift potential
 *  if flag == 0, do not shift
 *  if flag == 1, shift by potential and force value to zero at rCut
 */
void PairLJMulti::linearShift(const int flag    //!< initialization flag
  ) {
  cutShift(flag);
  linearShiftFlag_ = flag;
  peLinearShiftij_.resize(epsij_.size(), vector<double>(epsij_.size()));
  rCutij_.clear();
  rCutij_.resize(epsij_.size(), vector<double>(epsij_.size(), rCut_));
  if (flag == 0) {
    myFill(0., peLinearShiftij_);
  } else if (flag == 1) {
    for (unsigned int i = 0; i < epsij_.size(); ++i) {
      for (unsigned int j = 0; j < epsij_.size(); ++j) {
        linearShiftijset(i, j, flag);
      }
    }
  } else {
    ASSERT(0, "Unrecognized linearShift flag(" << flag << ").");
  }
}

/**
 * set i-j shifting
 */
void PairLJMulti::cutShiftijset(
  const int itype,
  const int jtype,
  const int flag) {
  cutShiftFlag_ = flag;
  peShiftij_.resize(epsij_.size(), vector<double>(epsij_.size()));
  if (flag == 0) {
    myFill(0., peShiftij_);
  } else if (flag == 1) {
    const double sig = sigij_[itype][jtype], rc = rCutij_[itype][jtype],
      eps = epsij_[itype][jtype];
    double rinv;
    if (sigrefFlag_ != 1) {
      rinv = sig/rc;
    } else {
      const double sigref = sigRefij_[itype][jtype];
      rinv = sigref / (rc - sig + sigref);
    }
    const double peShiftLJ = -4.*eps*(pow(rinv, 2*alpha_)
                                    - pow(rinv, alpha_));
    
    double peShiftY = 0;
    ASSERT(sigrefFlag_ == 0 || yukawa_ == 0, "sigref not implemented here");
    if (yukawa_ == 1) {
      peShiftY = -eps*yukawaA_*exp(-yukawaK_*rc/sig)/(rc/sig);
    } else if (yukawa_ == 2) {
      peShiftY = -yukawaA_*exp(-yukawaK_*rc/sig)/(rc/sig);
    }

    peShiftij_[itype][jtype] = peShiftLJ + peShiftY;
    peShiftij_[jtype][itype] = peShiftij_[itype][jtype];
    lrcFlag = 0;
  } else {
    ASSERT(0, "Unrecognized cutShift flag(" << flag << ").");
  }
}

/**
 * set i-j shifting
 */
void PairLJMulti::linearShiftijset(
  const int itype,
  const int jtype,
  const int flag) {
  cutShiftijset(itype, jtype, flag);
  linearShiftFlag_ = flag;
  peLinearShiftij_.resize(epsij_.size(), vector<double>(epsij_.size()));
  if (flag == 0) {
    myFill(0., peLinearShiftij_);
  } else if (flag == 1) {
    const double sig = sigij_[itype][jtype], rc = rCutij_[itype][jtype],
      eps = epsij_[itype][jtype];
    double rinv, sigtmp;
    if (sigrefFlag_ != 1) {
      rinv = sig/rc;
      sigtmp = sig;
    } else {
      const double sigref = sigRefij_[itype][jtype];
      rinv = sigref / (rc - sig + sigref);
      sigtmp = sigref;
    }
    
    const double peLShiftLJ = -4.*eps/sigtmp*(-2*alpha_*pow(rinv, (2*alpha_) 
      + 1) + alpha_*pow(rinv, alpha_ + 1));

    ASSERT(sigrefFlag_ == 0 || yukawa_ == 0, "sigref not implemented here");
    double peLShiftY = 0;
    if (yukawa_ == 1) {
      peLShiftY = eps*yukawaA_*exp(-yukawaK_*rc/sig)/(rc/sig)
        * (yukawaK_ + sig/rc) / sig;
    } else if (yukawa_ == 2) {
      peLShiftY = yukawaA_*exp(-yukawaK_*rc/sig)/(rc/sig)
        * (yukawaK_ + sig/rc) / sig;
    }

    peLinearShiftij_[itype][jtype] = peLShiftLJ + peLShiftY;
    peLinearShiftij_[jtype][itype] = peLinearShiftij_[itype][jtype];
  } else {
    ASSERT(0, "Unrecognized linearShift flag(" << flag << ").");
  }
}

/**
 * write restart file
 */
void PairLJMulti::writeRestart(const char* fileName) {
  PairLJ::writeRestart(fileName);
  std::ofstream file(fileName, std::ios_base::app);

  file << "# expType " << expType_ << endl;
  if (yukawa_ != 0) {
    file << "# yukawaFlag " << yukawa_ << endl;
    file << std::setprecision(std::numeric_limits<double>::digits10+2)
         << "# yukawaA " << yukawaA_ << endl;
    file << std::setprecision(std::numeric_limits<double>::digits10+2)
         << "# yukawaK " << yukawaK_ << endl;
  }

  for (unsigned int i = 0; i < peLinearShiftij_.size(); ++i) {
    for (unsigned int j = i; j < peLinearShiftij_.size(); ++j) {
      if (peLinearShiftij_[i][j] != 0) {
        file << std::setprecision(std::numeric_limits<double>::digits10+2)
             << "# peLinearShift" << i << j << " 1" << endl;
      } else if (peShiftij_[i][j] != 0) {
        file << std::setprecision(std::numeric_limits<double>::digits10+2)
             << "# peShift" << i << j << " 1" << endl;
      }
    }
  }

  if (gaussian_ == 1) {
    file << "# nGaussianParams " << gausParam_.size() << endl;
    for (unsigned int i = 0; i < gausParam_.size(); ++i) {
      file << std::setprecision(std::numeric_limits<double>::digits10+2)
           << "# gausParam" << i << "height " << gausParam_[i][0] << endl;
      file << std::setprecision(std::numeric_limits<double>::digits10+2)
           << "# gausParam" << i << "position " << gausParam_[i][1] << endl;
      file << std::setprecision(std::numeric_limits<double>::digits10+2)
           << "# gausParam" << i << "spread " << gausParam_[i][2] << endl;
    }
  }

  if (lambdaFlag_ == 1) {
    file << "# lambdaFlag " << lambdaFlag_ << endl;
    for (int i = 0; i < space_->nParticleTypes(); ++i) {
      for (int j = 0; j < space_->nParticleTypes(); ++j) {
        file << "# lambdaParami" << i << "j" << j << " " << lambda_[i][j]
             << endl;
      }
    }
  }
}

/**
 * initialize long-range correcitons
 */
void PairLJMulti::initLRC() {
  lrcFlag = 1;
  lrcPreCalc_.resize(epsij_.size(), vector<double>(epsij_.size()));
  ASSERT(expType_ != 1, "LRC not implemented for expType(" << expType_ << ")");
  ASSERT(sigrefFlag_ == 0, "sigref not implemented here");
  for (unsigned int i = 0; i < epsij_.size(); ++i) {
    for (unsigned int j = 0; j < epsij_.size(); ++j) {
      const double sig = sigij_[i][j];
      const double eps = epsij_[i][j];
      const double rc = rCutij_[i][j];
      if (expType_ == 0) {
        lrcPreCalc_[i][j] = (8./3.)*PI*eps*pow(sig, 3)*(pow(sig/rc, 9)/3.
          - pow(sig/rc, 3));
      } else {
        ASSERT(0, "Unrecognized expType(" << expType_ << ")");
      }
    }
  }
}

/**
 * set i-j to wca
 */
void PairLJMulti::initWCA(const int itype, const int jtype) {
  const double sig = sigij_[itype][jtype];
  double rc;
  if (sigrefFlag_ != 1) {
    rc = pow(2, 1./alpha_)*sig;
  } else {
    const double sigref = sigRefij_[itype][jtype];
    rc = pow(2, 1./alpha_)*sigref + sig - sigref;
  }
  rCutijset(itype, jtype, rc);
  cutShiftijset(itype, jtype, 1);
}

void PairLJMulti::multiPartEnerAtomCutInner(const double &r2, const int &itype,
  const int &jtype) {
  const double sigij = sigij_[itype][jtype];
  // cheap energy switches cut-off to sigmaij
  if ( (!cheapEnergy_) || (r2 < sigij*sigij) ) {
    double r2inv;
    double r = 0;
    double sigref = 0;
    double peLJ = 0;
    if (sigrefFlag_ != 1) {
      r2inv = sigij * sigij / r2;
    } else {
      sigref = sigRefij_[itype][jtype];
      r = sqrt(r2);
      r2inv = sigref/(r - sigij + sigref);
      r2inv = r2inv*r2inv;
    
      // inner hard sphere
      if (r <= sigij - sigref) {
        peLJ += std::numeric_limits<double>::max()/1e10;
        return;
      }
    }
    double r6inv = r2inv*r2inv*r2inv;
    if (expType_ == 1) {
      r6inv = r6inv*r6inv;
    } else if (expType_ == 2) {
      r6inv = pow(r2inv, 16.6755*0.5);
    } else if (expType_ == 3) {
      r6inv = pow(r2inv, 25);
    } else if (expType_ == 4) {
      r6inv = pow(r2inv, 64);
    } else if (expType_ == 5) {
      r6inv = pow(r2inv, 12);
    } else if (expType_ == 6) {
      r6inv = pow(r2inv, 9);
    }
    const double epsij = epsij_[itype][jtype];
    peLJ = epsij * (4. * (r6inv*(r6inv - 1.)) + peShiftij_[itype][jtype]);
    if (linearShiftFlag_) {
      if (sigrefFlag_ != 1) r = sqrt(r2);
      peLJ += peLinearShiftij_[itype][jtype] * (r - rCutij_[itype][jtype]);
    }

    if (lambdaFlag_ != 0) {
      const double lambda = lambda_[itype][jtype];
      double rwca;
      if (sigrefFlag_ != 1) {
        rwca = pow(2., 1./alpha_)*sigij;
      } else {
        rwca = pow(2., 1./alpha_)*sigref + sigij - sigref;
      }

      if (r2 < rwca*rwca) {
        peLJ += epsij*(1.-lambda) + (lambda-1.)*peShiftij_[itype][jtype];
        if (linearShiftFlag_) peLJ += (lambda-1.)*(r - rCutij_[itype][jtype])
                                      *peLinearShiftij_[itype][jtype];
      } else {
        peLJ *= lambda;
      }
      // cout << "peLJ " << peLJ << endl;
    }
    peSRone_ += peLJ;

    // yukawa
    if (yukawa_ == 1) {
      // if (!linearShiftFlag_) r = sqrt(r2);
      if (!linearShiftFlag_) r = sqrt(r2);
      peSRone_ += epsij*yukawaA_ * exp(-yukawaK_*r/sigij)/r*sigij;
    } else if (yukawa_ == 2) {
      if (!linearShiftFlag_) r = sqrt(r2);
      peSRone_ += yukawaA_ * exp(-yukawaK_*r/sigij)/r*sigij;
    }

    // gaussians
    if (gaussian_ == 1) {
      if ( (!linearShiftFlag_) && (!( (yukawa_ == 1) || (yukawa_ == 2) ) ) ) {
        r = sqrt(r2);
      }
      for (unsigned int ig = 0; ig < gausParam_.size(); ++ig) {
        peSRone_ += gausParam_[ig][0]*exp(-pow((r - gausParam_[ig][1])
          / gausParam_[ig][2], 2));
      }
    }
//    if (linearShiftFlag_) {
//      cout << "pepart " << epsij * (4. * (r6inv*(r6inv - 1.)) + peShiftij_[itype][jtype])+ static_cast<int>(linearShiftFlag_)*(peLinearShiftij_[itype][jtype] * (sqrt(r2) - rCutij_[itype][jtype])) + yukawaA_ * exp(-yukawaK_*sqrt(r2)/sigij)/sqrt(r2)*sigij << " r6inv " << r6inv << " eps " << epsij << " peshi " << peShiftij_[itype][jtype] << " lsf " << linearShiftFlag_ << " sig " << sigij_[itype][jtype] << " r2 " << r2 << " yuk " << yukawa_ << " yukaA " << yukawaA_ << " yukaK " << yukawaK_ << endl;
//    } else {
//      cout << "pepart " << epsij * (4. * (r6inv*(r6inv - 1.)) + peShiftij_[itype][jtype])+ yukawaA_ * exp(-yukawaK_*sqrt(r2)/sigij)/sqrt(r2)*sigij << " r6inv " << r6inv << " eps " << epsij << " peshi " << peShiftij_[itype][jtype] << " lsf " << linearShiftFlag_ << " sig " << sigij_[itype][jtype] << " r2 " << r2 << " yuk " << yukawa_ << " yukaA " << yukawaA_ << " yukaK " << yukawaK_ << endl;
//    }
  }
}

/**
 * inner loop for potential energy and forces of all particles
 */
void PairLJMulti::allPartEnerForceInner(const double &r2, const double &dx,
  const double &dy, const double &dz, const int &itype, const int &jtype,
  const int &iMol, const int &jMol) {
  const double sigij = sigij_[itype][jtype];
  const double r2inv = sigij*sigij/r2;
  double r6inv = r2inv*r2inv*r2inv;
  if (expType_ == 1) {
    r6inv = r6inv*r6inv;
  } else if (expType_ == 2) {
    r6inv = pow(r2inv, 16.6755*0.5);
  } else if (expType_ == 3) {
    r6inv = pow(r2inv, 25);
  } else if (expType_ == 4) {
    r6inv = pow(r2inv, 64);
  } else if (expType_ == 5) {
    r6inv = pow(r2inv, 12);
  } else if (expType_ == 6) {
    r6inv = pow(r2inv, 9);
  }
  const double epsij = epsij_[itype][jtype];
  multiPartEnerAtomCutInner(r2, itype, jtype);
  double fPart = 8.*alpha_*epsij*(r6inv*r2inv*(r6inv - 0.5)) / sigij / sigij;
  if (linearShiftFlag_) fPart -= peLinearShiftij_[itype][jtype]/sqrt(r2);
  fCOM_[iMol][0] += fPart * dx;
  fCOM_[iMol][1] += fPart * dy;
  fCOM_[jMol][0] -= fPart * dx;
  fCOM_[jMol][1] -= fPart * dy;
  if (dimen_ > 2) {
    fCOM_[iMol][2] += fPart * dz;
    fCOM_[jMol][2] -= fPart * dz;
  }
}

/** initialize exponential type
 *   type0 is 12-6
 *   type1 is 24-12
 *   type2 is 2alpha - alpha; alpha=16.6755
 *   type3 is 2alpha - alpha; alpha=50
 *   type4 is 2alpha - alpha; alpha=128
 *   type5 is 2alpha - alpha; alpha=24
 *   type6 is 2alpha - alpha; alpha=18
 */
void PairLJMulti::initExpType(const int type) {
  expType_ = type;
  if (expType_ == 0) {
    alpha_ = 6;
  } else if (expType_ == 1) {
    alpha_ = 12;
  } else if (expType_ == 2) {
    alpha_ = 16.6755;
  } else if (expType_ == 3) {
    alpha_ = 50;
  } else if (expType_ == 4) {
    alpha_ = 128;
  } else if (expType_ == 5) {
    alpha_ = 24;
  } else if (expType_ == 6) {
    alpha_ = 18;
  } else {
    ASSERT(0, "Unrecognized expType(" << expType_ << ")");
  }
}

/**
 * set the order parameter
 */
void PairLJMulti::setOrder(const double order) {
//  const double orderOld = order_;
  Pair::setOrder(order);

  if (orderName_.compare("MMsigma3") == 0) {
    initWCA(1, 2);
    initWCA(2, 2);

//    // scale the box such that V_ex/V ~ constant
//    //  V_ex is obtained approximately by fit to f(x)=a+b*order^c
//    const double a = 4.07988, b = 5.22454, c = 2.99222,
//                 vexold = a+b*pow(orderOld, c),
//                 vex =    a+b*pow(order_  , c);
//    if ( (orderOld > 0) && (vex/vexold > 0) ) {
//      space_->scaleDomain(vex/vexold);
//    }

  } else if (orderName_.compare("MMsigSheetVes") == 0) {
    initWCA(0, 1);
    initWCA(0, 2);
    initWCA(1, 1);
    initWCA(1, 2);
    initWCA(2, 2);
  } else if (orderName_.compare("lambda01") == 0) {
    setLambdaij(0, 1, order);
  } else if (orderName_.compare("angle0") != 0) {
    ASSERT(0, "Unrecognized orderName_(" << orderName_ << ") in setOrder");
  }
}

/**
 *  add a gaussian on the potential
 *   U(r) = height * exp ( - ( (r-position)/spread)^2 )
 */
void PairLJMulti::addGaussian(const double height, const double position,
  const double spread) {
  gaussian_ = 1;
  vector<double> param;
  param.push_back(height);
  param.push_back(position);
  param.push_back(spread);
  gausParam_.push_back(param);
}

/**
 * set lambda parameter
 */
void PairLJMulti::setLambdaij(const double iType, const double jType,
  const double lambda) {
  lambdaFlag_ = 1;
  lambda_.resize(epsij_.size(), vector<double>(epsij_.size()));
  lambda_[iType][jType] = lambda;
  lambda_[jType][iType] = lambda;
}

/**
 * return the lrc contribution of one particle
 */
double PairLJMulti::computeLRC(const int ipart) {
  double enlrc = 0;
  const int iType = space_->type()[ipart];

  // search all particle types, and sum lrcs of each
  for (int jType = 0; jType < space_->nParticleTypes(); ++jType) {
    int n = space_->nType()[jType];
    enlrc += static_cast<double>(n)/space_->vol() * lrcPreCalc_[iType][jType];
  }
  return enlrc;
}

}  // namespace feasst

