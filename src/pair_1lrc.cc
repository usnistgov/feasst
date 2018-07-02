/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./pair_1lrc.h"

namespace feasst {

PairLRC::PairLRC(Space* space, const argtype &args)
  : Pair(space, args) {
  defaultConstruction_();
}

PairLRC::PairLRC(Space* space, const char* fileName)
  : Pair(space, fileName) {
  defaultConstruction_();
  string strtmp = fstos("linearShiftFlag", fileName);
  if (!strtmp.empty()) {
    linearShift(stoi(strtmp));
  } else {
    string strtmp = fstos("cutShiftFlag", fileName);
    if (!strtmp.empty()) cutShift(stoi(strtmp));
  }

  // set exponential type
  string str = fstos("expType", fileName);
  if (!str.empty()) {
    const int expType = stoi(str);
    if (expType == -1) {
      str = fstos("alphaparam", fileName);
      ASSERT(!str.empty(), "alpha parameter must be provided if expType==-1");
      initAlpha(stod(str));
    } else {
      initExpType(stoi(str));
    }
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

  // "global rcut" check first
  str = fstos("cutShiftFlag", fileName);
  if (!str.empty()) cutShift(stoi(str));
  if (cutShiftFlag_ == 0) cutShift(0);

  // overwrite global rcut with ij pairs if specified
  for (int i = 0; i < space_->nParticleTypes(); ++i) {
    for (int j = i; j < space_->nParticleTypes(); ++j) {  // only do j >= i
      stringstream ss;
      ss << "rCut" << i << "j" << j;
      string strtmp = fstos(ss.str().c_str(), fileName);
      if (!strtmp.empty()) {
        double rc_ij = fstod(ss.str().c_str(), fileName);
        rCutijset(i, j, rc_ij);
      }
    }
  }

  // once r_cuts are set in the correct order, shift/cut as necessary
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

  // lrc flag is read in Pair
  if (lrcFlag == 1) {
    initLRC();
  }
}

void PairLRC::defaultConstruction_() {
  className_.assign("PairLRC");
  linearShift(0);
  initExpType(0);
  yukawa_ = 0;
  yukawaA_ = 0;
  yukawaK_ = 0;
}

void PairLRC::initLRC() {
  lrcFlag = 1;
  lrcPreCalc_.resize(epsij_.size(), vector<double>(epsij_.size()));
  ASSERT(sigrefFlag_ == 0, "sigref not implemented here");
  for (unsigned int i = 0; i < epsij_.size(); ++i) {
    for (unsigned int j = 0; j < epsij_.size(); ++j) {
      const double sig = sigij_[i][j];
      const double eps = epsij_[i][j];
      double rc = 0.;
      if (rCutij_.size() == 0 || atomCut_ == 0) {
        rc = rCut_;
      } else {
        rc = rCutij_[i][j];
        // cout << "rc" << i << j << " " << rc << endl;
      }
      if ( (fabs(sig) > DTOL) && (fabs(eps) > DTOL) ) {
        lrcPreCalc_[i][j] = (8./3.)*PI*eps*pow(sig, 3)*(pow(sig/rc, 9)/3.
          - pow(sig/rc, 3));
        lrcPreCalc_[j][i] = lrcPreCalc_[i][j];
      }
    }
  }
//  cout << "Init LRCs" << endl << vec2str(lrcPreCalc_) << endl;
}

double PairLRC::computeLRC(const vector<int> msite) {
  double enlrc = 0.;
  if (lrcFlag == 0) return enlrc;
  if (lrcPreCalc_.size() == 0) {
    initLRC();
  }
  // If msite is empty, loop over all sites
  if (msite.size() == 0) {
    for (int iType = 0; iType < space_->nParticleTypes(); ++iType) {
      for (int jType = 0; jType < space_->nParticleTypes(); ++jType) {
        const int nanb = space_->nType()[iType]*space_->nType()[jType];
        enlrc += static_cast<double>(nanb)/space_->volume()
                *lrcPreCalc_[iType][jType];
      }
    }
  } else {
    // Find the number of particles of each type
    vector<int> msiteTypes(space_->nParticleTypes(), 0);
    for (int imsite = 0; imsite < static_cast<int>(msite.size()); ++imsite) {
      ++msiteTypes[space()->type()[msite[imsite]]];
    }

    for (int iType = 0; iType < space_->nParticleTypes(); ++iType) {
      for (int jType = 0; jType < space_->nParticleTypes(); ++jType) {
        const double ni = space_->nType()[iType],
                     nia = msiteTypes[iType],
                     nj = space_->nType()[jType],
                     nja = msiteTypes[jType];
        enlrc += (ni*nja + nia*nj - nia*nja)/space_->volume()
                 *lrcPreCalc_[iType][jType];
      }
    }
  }
  return enlrc;
}

void PairLRC::linearShift(const int flag) {
  linearShiftFlag_ = flag;
  cutShift(flag);
  peLinearShiftij_.resize(epsij_.size(), vector<double>(epsij_.size()));
  rCutij_.clear();
  rCutij_.resize(epsij_.size(), vector<double>(epsij_.size(), rCut_));
  if (flag == 0) {
    fill(0., peLinearShiftij_);
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

void PairLRC::cutShift(const int flag) {
  cutShiftFlag_ = flag;
  peShiftij_.resize(epsij_.size(), vector<double>(epsij_.size()));
  rCutij_.clear();
  rCutij_.resize(epsij_.size(), vector<double>(epsij_.size(), rCut_));
  fill(rCut_, rCutij_);
  if (flag == 0) {
    fill(0., peShiftij_);
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

void PairLRC::cutShiftijset(
  const int itype,
  const int jtype,
  const int flag) {
  cutShiftFlag_ = flag;
  peShiftij_.resize(epsij_.size(), vector<double>(epsij_.size()));
  if (flag == 0) {
    fill(0., peShiftij_);
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
    if (yukawa_ == 1) {
      peShiftY = -eps*yukawaA_*exp(-yukawaK_*rc/sig)/(rc/sig);
    } else if (yukawa_ == 2) {
      peShiftY = -yukawaA_*exp(-yukawaK_*rc/sig)/(rc/sig);
    } else if (yukawa_ == 3) {
      peShiftY = -yukawaAij_[itype][jtype]
                 *exp(-yukawaKij_[itype][jtype]*rc/sig)/(rc/sig);
    }

    peShiftij_[itype][jtype] = peShiftLJ + peShiftY;
    peShiftij_[jtype][itype] = peShiftij_[itype][jtype];
    lrcFlag = 0;
  } else {
    ASSERT(0, "Unrecognized cutShift flag(" << flag << ").");
  }
}

void PairLRC::linearShiftijset(
  const int itype,
  const int jtype,
  const int flag) {
  cutShiftijset(itype, jtype, flag);
  linearShiftFlag_ = flag;
  peLinearShiftij_.resize(epsij_.size(), vector<double>(epsij_.size()));
  if (flag == 0) {
    fill(0., peLinearShiftij_);
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

    double peLShiftY = 0;
    if (yukawa_ == 1) {
      peLShiftY = eps*yukawaA_*exp(-yukawaK_*rc/sig)/(rc/sig)
        * (yukawaK_ + sig/rc) / sig;
    } else if (yukawa_ == 2) {
      peLShiftY = yukawaA_*exp(-yukawaK_*rc/sig)/(rc/sig)
        * (yukawaK_ + sig/rc) / sig;
    } else if (yukawa_ == 3) {
      peLShiftY = yukawaAij_[itype][jtype]
        * exp(-yukawaKij_[itype][jtype]*rc/sig)/(rc/sig)
        * (yukawaKij_[itype][jtype] + sig/rc) / sig;
    }

    peLinearShiftij_[itype][jtype] = peLShiftLJ + peLShiftY;
    peLinearShiftij_[jtype][itype] = peLinearShiftij_[itype][jtype];
  } else {
    ASSERT(0, "Unrecognized linearShift flag(" << flag << ").");
  }
}

void PairLRC::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  if (linearShiftFlag_ != 0) {
    file << "# linearShiftFlag " << linearShiftFlag_ << endl;
  } else if (cutShiftFlag_ != 0) {
    file << "# cutShiftFlag " << cutShiftFlag_ << endl;
  }

  file << "# expType " << expType_ << endl;
  file << "# alphaparam " << alpha_ << endl;
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

  if (lrcFlag == 1) {
    file << "# lrcFlag " << lrcFlag << endl;
  }
}

void PairLRC::initExpType(const int type) {
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
  } else if (expType_ == -1) {
    ASSERT(0, "Initialize expType=-1 via initAlpha instead");
  } else {
    ASSERT(0, "Unrecognized expType(" << expType_ << ")");
  }
}

void PairLRC::initAlpha(const double alpha) {
  alpha_ = alpha;
  expType_ = -1;
}

}  // namespace feasst
