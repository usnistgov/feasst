/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./pair_lj.h"
#include "./arguments.h"

namespace feasst {

PairLJ::PairLJ(Space* space, const argtype &args)
  : PairLRC(space, args) {
  defaultConstruction_();
  argparse_.initArgs(className_, args);

  // parse molType
  std::stringstream ss;
  ss << install_dir() << "/forcefield/data.lj";
  std::string molType(ss.str());
  if (!argparse_.key("molType").empty()) {
    molType = argparse_.str();
  // if molType is not provided, try optional forcefield arg
  } else if (!argparse_.key("molTypeInForcefield").empty()) {
    ss.str("");
    ss << install_dir() << "/forcefield/" << argparse_.str();
    molType = ss.str();
  }
  if (molType != "none") {
    initData(molType);

    // apply cutType
    const std::string cutType = argparse_.key("cutType").dflt("lrc").str();
    if (cutType == "lrc") {
      cutShift(0);
      initLRC();
    } else if (cutType == "cutShift") {
      cutShift(1);
    } else if (cutType == "linearShift") {
      linearShift(1);
    } else if (cutType == "none") {
      lrcFlag = 0;
      cutShift(0);
    } else {
      ASSERT(0, "unrecognized cutType(" << cutType << ")");
    }

    // initialize potential energy
    initEnergy();
  }
  argparse_.checkAllArgsUsed();
}

PairLJ::PairLJ(Space* space, const char* fileName)
  : PairLRC(space, fileName) {
  defaultConstruction_();

  std::string str = fstos("nGaussianParams", fileName);
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
}

void PairLJ::defaultConstruction_() {
  className_.assign("PairLJ");
  lambdaFlag_ = 0;
  gaussian_ = 0;
}

void PairLJ::initEnergy() {
  // zero accumulators: potential energy, force, and virial
  std::fill(pe_.begin(), pe_.end(), 0.);
  fill(0., f_);
  fill(0., vr_);
  peLJ_ = 0;
  fCOM_.clear();
  fCOM_.resize(space_->nMol(), vector<double>(dimen_, 0.));
  peTot_ = allPartEnerForce(2);
  peLJ_ = peSRone_;
  peLRC_ = peLRCone_;
}

double PairLJ::allPartEnerForce(const int flag) {
  peSRone_ = 0;
  // standard long range corrections
  peLRCone_ = computeLRC();

  if (flag == 0) {
    peSRone_ = peTot() - peLRCone_;
    return peSRone_ + peLRCone_;
  } else {
    // zero accumulators: potential energy and force
    fCOM_.clear();
    fCOM_.resize(space_->nMol(), vector<double>(dimen_, 0.));

    int noCell = 0;
    if (flag == 2) noCell = 1;
    return pairLoopSite_(noCell) + peLRCone_;
  }
  ASSERT(0, "cell type(" << space_->cellType() << ")");
  return 1e300;
}

double PairLJ::multiPartEner(const vector<int> mpart, const int flag) {
  if (flag == 0) {}  // remove unused parameter warning

  // zero potential energy contribution of particle ipart
  peSRone_ = 0;
  peLRCone_ = 0;
  if (!cheapEnergy_) {
    peLRCone_ += computeLRC(mpart);
  }
  return pairLoopSite_(mpart) + peLRCone_;
}

void PairLJ::writeRestart(const char* fileName) {
  PairLRC::writeRestart(fileName);
  std::ofstream file(fileName, std::ios_base::app);
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

void PairLJ::initLRC() {
  ASSERT(expType_ == 0, "LRC not implemented for expType(" << expType_ << ")");
  PairLRC::initLRC();
}

void PairLJ::initWCA(const int itype, const int jtype) {
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

void PairLJ::initWCA() {
  for (unsigned int i = 0; i < epsij_.size(); ++i) {
    for (unsigned int j = 0; j < epsij_.size(); ++j) {
      initWCA(i, j);
    }
  }
}

void PairLJ::setOrder(const double order) {
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

void PairLJ::addGaussian(const double height, const double position,
  const double spread) {
  gaussian_ = 1;
  vector<double> param;
  param.push_back(height);
  param.push_back(position);
  param.push_back(spread);
  gausParam_.push_back(param);
}

void PairLJ::setLambdaij(const double iType, const double jType,
  const double lambda) {
  lambdaFlag_ = 1;
  lambda_.resize(epsij_.size(), vector<double>(epsij_.size()));
  lambda_[iType][jType] = lambda;
  lambda_[jType][iType] = lambda;
}

void PairLJ::pairSiteSite_(const int &iSiteType, const int &jSiteType,
  double * energy, double * force, int * neighbor, const double &dx,
  const double &dy, const double &dz) {
  const double r2 = dx*dx + dy*dy + dz*dz;
  *energy = 0;
  *neighbor = 0;
  *force = 0;
  const double sigij = sigij_[iSiteType][jSiteType];
  // cheap energy switches cut-off to sigmaij
  if ( (!cheapEnergy_) || (r2 < sigij*sigij) ) {
    double r2inv;
    double r = 0;
    double sigref = 0;
    double peLJ = 0;
    if (sigrefFlag_ != 1) {
      r2inv = sigij * sigij / r2;
    } else {
      sigref = sigRefij_[iSiteType][jSiteType];
      r = sqrt(r2);
      r2inv = sigref/(r - sigij + sigref);
      r2inv = r2inv*r2inv;

      // inner hard sphere
      if (r <= sigij - sigref) {
        peLJ += NUM_INF;
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
    } else if (expType_ == -1) {
      r6inv = pow(r2inv, 0.5*alpha_);
    }
    const double epsij = epsij_[iSiteType][jSiteType];
    peLJ += epsij * 4. * (r6inv*(r6inv - 1.)) + peShiftij_[iSiteType][jSiteType];
    if (linearShiftFlag_) {
      if (sigrefFlag_ != 1) r = sqrt(r2);
      peLJ += peLinearShiftij_[iSiteType][jSiteType] * (r - rCutij_[iSiteType][jSiteType]);
    }

    *force = 8.*alpha_*(r6inv*r2inv*(r6inv - 0.5));
    if (linearShiftFlag_) {
      *force -= peLinearShiftij_[iSiteType][jSiteType]/sqrt(r2);
    }

    if (lambdaFlag_ != 0) {
      const double lambda = lambda_[iSiteType][jSiteType];
      double rwca;
      if (sigrefFlag_ != 1) {
        rwca = pow(2., 1./alpha_)*sigij;
      } else {
        rwca = pow(2., 1./alpha_)*sigref + sigij - sigref;
      }

      if (r2 < rwca*rwca) {
        peLJ += epsij*(1.-lambda) + (lambda-1.)*peShiftij_[iSiteType][jSiteType];
        if (linearShiftFlag_) peLJ += (lambda-1.)*(r - rCutij_[iSiteType][jSiteType])
                                      *peLinearShiftij_[iSiteType][jSiteType];
      } else {
        peLJ *= lambda;
      }
      // cout << "peLJ " << peLJ << endl;
    }
    *energy = peLJ;

    // yukawa
    if (yukawa_ == 1) {
      // if (!linearShiftFlag_) r = sqrt(r2);
      if (!linearShiftFlag_) r = sqrt(r2);
      *energy += epsij*yukawaA_ * exp(-yukawaK_*r/sigij)/r*sigij;
    } else if (yukawa_ == 2) {
      if (!linearShiftFlag_) r = sqrt(r2);
      *energy += yukawaA_ * exp(-yukawaK_*r/sigij)/r*sigij;
    } else if (yukawa_ == 3) {
      if (!linearShiftFlag_) r = sqrt(r2);
      *energy += yukawaAij_[iSiteType][jSiteType]
                 * exp(-yukawaKij_[iSiteType][jSiteType]*r/sigij)/r*sigij;
    }

    // gaussians
    if (gaussian_ == 1) {
      if ( (!linearShiftFlag_) && (!( (yukawa_ == 1) || (yukawa_ == 2) ) ) ) {
        r = sqrt(r2);
      }
      for (unsigned int ig = 0; ig < gausParam_.size(); ++ig) {
        *energy += gausParam_[ig][0]*exp(-pow((r - gausParam_[ig][1])
          / gausParam_[ig][2], 2));
      }
    }
//    if (linearShiftFlag_) {
//      cout << "pepart " << epsij * (4. * (r6inv*(r6inv - 1.)) + peShiftij_[iSiteType][jSiteType])+ static_cast<int>(linearShiftFlag_)*(peLinearShiftij_[iSiteType][jSiteType] * (sqrt(r2) - rCutij_[iSiteType][jSiteType])) << " r6inv " << r6inv << " eps " << epsij << " peshi " << peShiftij_[iSiteType][jSiteType] << " lsf " << linearShiftFlag_ << " sig " << sigij_[iSiteType][jSiteType] << " r2 " << r2 << " yuk " << yukawa_ << " yukaA " << yukawaA_ << " yukaK " << yukawaK_ << endl;
//    //+ yukawaA_ * exp(-yukawaK_*sqrt(r2)/sigij)/sqrt(r2)*sigij << " r6inv " << r6inv << " eps " << epsij << " peshi " << peShiftij_[iSiteType][jSiteType] << " lsf " << linearShiftFlag_ << " sig " << sigij_[iSiteType][jSiteType] << " r2 " << r2 << " yuk " << yukawa_ << " yukaA " << yukawaA_ << " yukaK " << yukawaK_ << endl;
//  }
//    } else {
//      cout << "pepart " << epsij * (4. * (r6inv*(r6inv - 1.)) + peShiftij_[iSiteType][jSiteType])+ yukawaA_ * exp(-yukawaK_*sqrt(r2)/sigij)/sqrt(r2)*sigij << " r6inv " << r6inv << " eps " << epsij << " peshi " << peShiftij_[iSiteType][jSiteType] << " lsf " << linearShiftFlag_ << " sig " << sigij_[iSiteType][jSiteType] << " r2 " << r2 << " yuk " << yukawa_ << " yukaA " << yukawaA_ << " yukaK " << yukawaK_ << endl;
//    }
    *neighbor = 1;
  }
}

double PairLJ::pairLoopSite_(
  const vector<int> &siteList,
  const int noCell) {
  if ( (siteList.size() != 1) ||
       (epsij_.size() > 1) ||
       (space_->tilted()) ||
       (dimen_ != 3) ||
       (yukawa_ != 0) ||
       (expType_ != 0) ||
       (lambdaFlag_ != 0) ||
       (gaussian_ != 0) ) {
    return Pair::pairLoopSite_(siteList, noCell);
  }

  // shorthand for read-only space variables
  const int natom = space_->natom();
  const vector<double> &x = space_->x();
  const vector<int> &mol = space_->mol();
  const vector<double> &boxLength = space_->boxLength();

  // declare variables for optimization
  double r6inv, r2, xi, yi, zi, dx, dy, dz;
  const double rCut = rCutij_[0][0];
  const double rCutSq = rCut*rCut;
  const double peShift = peShiftij_[0][0];
  const double epsilon = epsij_[0][0];
  const double sigma = sigij_[0][0];
  const double sigmaSq = sigma*sigma;
  double peLinearShift = 0.;
  if (linearShiftFlag_) {
    peLinearShift = peLinearShiftij_[0][0];
  }

  // PBC optimization variables
  const double lx = boxLength[0];
  const double ly = boxLength[1];
  const double lz = boxLength[2];
  const double halflx = lx/2., halfly = ly/2., halflz = lz/2.;
  const int xpbc = space_->periodic(0),
            ypbc = space_->periodic(1),
            zpbc = space_->periodic(2);

  initNeighCutPEMap(siteList);

  // loop through all particles in siteList
  for (unsigned int ii = 0; ii < siteList.size(); ++ii) {
    const int ipart = siteList[ii];
    const int iMol = mol[ipart];
    xi = x[dimen_*ipart];
    yi = x[dimen_*ipart+1];
    zi = x[dimen_*ipart+2];

    // loop through all particles interacting with ipart
    for (int jpart = 0; jpart < natom; ++jpart) {
      if (iMol != mol[jpart]) {
        // separation distance with periodic boundary conditions
        dx = xi - x[dimen_*jpart];
        dy = yi - x[dimen_*jpart+1];
        dz = zi - x[dimen_*jpart+2];
        if (xpbc && dx >  halflx) dx -= lx;
        if (xpbc && dx < -halflx) dx += lx;
        if (ypbc && dy >  halfly) dy -= ly;
        if (ypbc && dy < -halfly) dy += ly;
        if (zpbc && dz >  halflz) dz -= lz;
        if (zpbc && dz < -halflz) dz += lz;
        r2 = dx*dx + dy*dy + dz*dz;

        // no interaction beyond cut-off distance
        if (r2 < rCutSq) {
          const double r2inv = sigmaSq/r2;
          r6inv = r2inv*r2inv*r2inv;
          peSRone_ += 4. * epsilon* (r6inv*(r6inv - 1.)) + peShift;
          if (linearShiftFlag_) {
              peSRone_ += peLinearShift * (sqrt(r2) - rCut);
          }
          setNeighbor_(r2, ii, jpart, 0, 0);
        }
      }
    }
  }
  //return peSRone_;
  if (peMapOn_ == 1) {
    peSRone_ = peSRoneAlt_;
  }
  return peSRone_;
}

void PairLJ::update(
  const vector<int> mpart,
  const int flag,
  const char* uptype) {
  if (neighOn_) {
    // rebuilt neighlist if all particles are updated
    if (static_cast<int>(mpart.size()) != space_->natom()) {
      updateBase(mpart, flag, uptype, neigh_, neighOne_, neighOneOld_);
  //  updateBase(mpart, flag, uptype, neighCut_, neighCutOne_, neighCutOneOld_);
    }
  }
  std::string uptypestr(uptype);

  if (uptypestr.compare("store") == 0) {
    if (flag == 0 || flag == 2 || flag == 3) {
      deLJ_ = peSRone_;
      deLRC_ = peLRCone_;
//      cout << "storing deLJ_ " << deLJ_ << " deLRC_ " << deLRC_ << endl;
    }
  }

  if (uptypestr.compare("update") == 0) {
    if (flag == 0) {
      peLJ_ += peSRone_ - deLJ_;
      peLRC_ += peLRCone_ - deLRC_;
      peTot_ += peSRone_ - deLJ_ + peLRCone_ - deLRC_;
//      cout << "de " << peSRone_ - deLJ_ << endl;
//      cout << "delrc " << peLRCone_ - deLRC_ << endl;
    }
    if (flag == 2) {
      peLJ_ -= deLJ_;
      peLRC_ -= deLRC_;
      peTot_ -= deLJ_ + deLRC_;
    }
    if (flag == 3) {
      peLJ_ += deLJ_;
      peLRC_ += deLRC_;
      peTot_ += deLJ_ + deLRC_;
    }
  }
}

void PairLJ::initScreenedElectroIJ(const int itype, const int jtype,
  const double A, const double K) {
  yukawa_ = 3;
  if (static_cast<int>(yukawaAij_.size()) <= itype) yukawaAij_.resize(itype+1);
  if (static_cast<int>(yukawaAij_.size()) <= jtype) yukawaAij_.resize(jtype+1);
  if (static_cast<int>(yukawaKij_.size()) <= itype) yukawaKij_.resize(itype+1);
  if (static_cast<int>(yukawaKij_.size()) <= jtype) yukawaKij_.resize(jtype+1);
  if (static_cast<int>(yukawaAij_[itype].size()) <= jtype) {
    yukawaAij_[itype].resize(jtype+1);
  }
  if (static_cast<int>(yukawaAij_[jtype].size()) <= itype) {
    yukawaAij_[jtype].resize(itype+1);
  }
  if (static_cast<int>(yukawaKij_[itype].size()) <= jtype) {
    yukawaKij_[itype].resize(jtype+1);
  }
  if (static_cast<int>(yukawaKij_[jtype].size()) <= itype) {
    yukawaKij_[jtype].resize(itype+1);
  }
  yukawaAij_[itype][jtype] = A;
  yukawaAij_[jtype][itype] = A;
  yukawaKij_[itype][jtype] = K;
  yukawaKij_[jtype][itype] = K;
}

shared_ptr<PairLJ> makePairLJ(Space* space, const argtype &args) {
  return make_shared<PairLJ>(space, args);
}

shared_ptr<PairLJ> makePairLJ(shared_ptr<Space> space, const argtype &args) {
  return makePairLJ(space.get(), args);
}

}  // namespace feasst
