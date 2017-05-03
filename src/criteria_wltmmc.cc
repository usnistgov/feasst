/**
 * \file
 *
 * \brief acceptance criteria for Wang-Landau monte carlo trials
 */

#include "./criteria_wltmmc.h"
#include "./space.h"
#include "./pair.h"
#include "./mins.h"

namespace feasst {

/**
 * Constructor
 */
CriteriaWLTMMC::CriteriaWLTMMC(const double beta,  //!< inverse temperature
  const double activ,  //!< activity
  const char* mType,  //!< definition of macrostate
  const double mMin,  //!< minimum value of macrostate
  const double mMax,  //!< maximum value of macrostate
  const int nBin)     //!< number of bins for macrostate
  : Criteria(beta, activ) {
  defaultConstruction();
  mType_ = mType;
  initBins(mMax, mMin, nBin);
  zeroStat();
}

/**
 * construct from file
 */
CriteriaWLTMMC::CriteriaWLTMMC(const char* fileName
  )
    : Criteria(fileName) {
  defaultConstruction();
  mType_.assign(fstos("mType", fileName));
  mMin_ = fstod("mMin", fileName);
  mMax_ = fstod("mMax", fileName);
  nBin_ = fstoi("nBin", fileName);
  initBins(mMax_, mMin_, nBin_);
  wlFlat_ = fstoi("wlFlat", fileName);
  lnf_ = fstod("lnf", fileName);
  g_ = fstod("gwlmod", fileName);
  nSweep_ = fstoi("nSweep", fileName);
  collect_ = fstoi("collect", fileName);
  tmmc_ = fstoi("tmmc", fileName);
  lnfCollect_ = fstod("lnfCollect", fileName);
  lnfTMMC_ = fstod("lnfTMMC", fileName);
  wlFlatFactor_ = fstod("wlFlatFactor", fileName);
  nSweepVisPerBin_ = fstoi("nSweepVisPerBin", fileName);
  readlnPIEnerCol(fileName);
}

/**
 * defaults in constructor
 */
void CriteriaWLTMMC::defaultConstruction() {
  className_.assign("CriteriaWLTMMC");
  g_ = 0.5;
  lnfCollect_ = 1e-6;
  lnfTMMC_ = 1e-99;
  wlFlatFactor_ = 0.8;
  nSweepVisPerBin_ = 100;
  collect_ = false;
  tmmc_ = false;
  phaseBoundary_ = 0;
  nSmooth_ = 10;
}

/**
 * clone design pattern
 */
CriteriaWLTMMC* CriteriaWLTMMC::clone() const {
  CriteriaWLTMMC* c = new CriteriaWLTMMC(*this);
  c->reconstruct();
  return c;
}
shared_ptr<CriteriaWLTMMC> CriteriaWLTMMC::cloneShrPtr() const {
  return(std::static_pointer_cast<CriteriaWLTMMC, Criteria>(cloneImpl_()));
}
shared_ptr<Criteria> CriteriaWLTMMC::cloneImpl_() const {
  shared_ptr<CriteriaWLTMMC> c = make_shared<CriteriaWLTMMC>(*this);
  c->reconstruct();
  return c;
}

/**
 * acceptance criteria for trial moves
 *  returns 1 if accepted, 0 otherwise
 */
int CriteriaWLTMMC::accept(
  const double lnpMet,   //!< log of the metropolis acceptance probability
  const double peNew,    //!< potential energy of proposed configuration
  const char* moveType,  //!< type of move
  const int reject) {    //!< outright reject move if 1
  int returnVal = -1;
  mNew_ = mMin_ - 1;
  if ( (mType_.compare("nmol") == 0) ||
       (mType_.compare("nmolstage") == 0) ) {
    std::string moveTypeStr(moveType);
    if (moveTypeStr.compare("move") == 0) {
      mNew_ = mOld_;
    } else if (moveTypeStr.compare("add") == 0) {
      mNew_ = mOld_ + mBin_;
    } else if (moveTypeStr.compare("del") == 0) {
      mNew_ = mOld_ - mBin_;
    } else {
      ASSERT(0, "move type (" << moveTypeStr
        << ") is not recognized in wltmmc acceptance criteria");
    }
  } else if (mType_.compare("energy") == 0) {
    mNew_ = peNew;
  } else if ( (mType_.compare("pairOrder") == 0) ||
              (mType_.compare("pressure") == 0) ||
              (mType_.compare("lnpres") == 0) ||
              (mType_.compare("beta") == 0) ) {
    std::string moveTypeStr(moveType);
    if ( (moveTypeStr.compare("move") == 0) ||
         (moveTypeStr.compare("add")  == 0) ||
         (moveTypeStr.compare("del")  == 0) ) {
      mNew_ = mOld_;
    } else {
      mNew_ = std::stod(moveType);
    }
  } else {
    ASSERT(0, "unrecognized macrostate type (" << mType_ << ")");
  }
  const int mOldBin = bin(mOld_);
  mout_("verbose", std::ostringstream().flush() << "mNew " << mNew_);

  if (cTripleBanded_) {
    double pMet = exp(lnpMet);
    const int mNewBin = bin(mNew_);
    if ( (mNew_ > mMax_) || (mNew_ < mMin_) || (reject == 1) ) {
      returnVal = 0;
      if (reject == 1) pMet = 0;
    } else if (uniformRanNum() < exp( lnPI_[mOldBin]
                                    - lnPI_[mNewBin] + lnpMet)) {
      returnVal = 1;
    } else {
      returnVal = 0;
    }
  
    // update pe_ vector
    if (collect_) {
      if (returnVal == 0) {
        pe_[mOldBin].accumulate(peOld_);
      } else {
        pe_[mNewBin].accumulate(peNew);
      }
    }

    mUpdate(mOldBin, mNewBin, pMet, returnVal);
  }
  if (returnVal == 0) mNew_ = mOld_;
  ASSERT( (returnVal == 0) || (returnVal == 1), "returnVal(" << returnVal
    << "not 0 or 1");
  return returnVal;
}

/**
 * store macrostate variables of old configuration
 */
void CriteriaWLTMMC::store(const Space* space, Pair* pair) {
  if (mType_.compare("nmol") == 0) {
    mOld_ = space->nMol();
  } else if (mType_.compare("nmolstage") == 0) {
    mOld_ = space->nMol();
    if (space->tagStage() != 0) mOld_ += space->tagStage() - 1;
  } else if (mType_.compare("energy") == 0) {
    mOld_ = pair->peTot();
  } else if (mType_.compare("pairOrder") == 0) {
    mOld_ = pair->order();
  } else if (mType_.compare("beta") == 0) {
    mOld_ = beta_;
  } else if (mType_.compare("pressure") == 0) {
    mOld_ = pressure_;
  } else if (mType_.compare("lnpres") == 0) {
    mOld_ = log(pressure_);
  }
  peOld_ = pair->peTot();
  mout_("verbose", std::ostringstream().flush() << "mold " << mOld_);
  ASSERT((mOld_ <= mMax_) && (mOld_ >= mMin_), "current macrostate variable ("
    << mOld_ << ") is beyond the limits (" << mMin_ << " to " << mMax_ << ")");
  if (tmmc_ == false) flatCheck();
}

/**
 * if flatness criteria is met, reset histogram and reduce update factor
 */
int CriteriaWLTMMC::flatCheck() {
  // flatness criteria is met when the minimum value of the histogram is within
  // a factor, flatFac, of the average histogram value
  if (*std::min_element(h_.begin(), h_.end()) > wlFlatFactor_ * myVecAv(h_)) {
    std::fill(h_.begin(), h_.end(), 0);
    lnf_ *= g_;
    ++wlFlat_;
    if ( (lnf_ < lnfCollect_) && (!collect_) ) collectInit();
    if ( (lnf_ < lnfTMMC_) && (!tmmc_) ) tmmcInit();
    return 1;
  }
  return 0;
}

/**
 * update WL and/or the collection matrix
 */
void CriteriaWLTMMC::mUpdate(const int mOldBin,   //!< bin of old state
  const int mNewBin,   //!< bin of new state
  const double pmet,  //!< transition probability
  const int acceptFlag  //!< transition accepted if == 1
  ) {
  if (collect_) {
    const double p = std::min(1., pmet);
    if (cTripleBanded_) {
      int index = -1;
      if (mOldBin == mNewBin) {
        index = 1;
      } else {
        if (mOldBin - 1 == mNewBin) {
          index = 0;
        } else if (mOldBin + 1 == mNewBin) {
          index = 2;
        } else {
          ASSERT(0, "mOldBin = " << mOldBin << " while mNewBin = " << mNewBin);
        }

        // keep track of visited states for sweeps
        if ( (tmmc_ == true) && (acceptFlag == 1) ) {
          ++h_[mNewBin];
        }
      }
      C_[mOldBin][index] += p;
      C_[mOldBin][1] += 1-p;
    } else {
      C_[mOldBin][mNewBin] += p;
      C_[mOldBin][mOldBin] += 1-p;
    }
  }

  // Wang-Landau update
  if (tmmc_ == false) {
    int bin = mOldBin;
    if (acceptFlag == 1) bin = mNewBin;
    ++h_[bin];
    lnPI_[bin] += lnf_;
  }

  // number of tunnels
  if (mNewBin == mMin_) {
    if (nTunnelPrev_ == 1) {
      ++nTunnels_;
      nTunnelPrev_ = 0;
    } else if (nTunnelPrev_ == -1) {
      nTunnelPrev_ = 0;
    }
  } else if (mNewBin == mMax_) {
    if (nTunnelPrev_ == 0) {
      ++nTunnels_;
      nTunnelPrev_ = 1;
    } else if (nTunnelPrev_ == -1) {
      nTunnelPrev_ = 1;
    }
  }
}

/**
 * update the collection matrix
 */
void CriteriaWLTMMC::lnPIupdate() {
  if (tmmc_) {
    c2lnPI(C_, &lnPI_);

    // update number of sweeps
    if (*std::min_element(h_.begin(), h_.end()) >= nSweepVisPerBin_) {
      ++nSweep_;
      std::fill(h_.begin(), h_.end(), 0);
    }
  }
}

/**
 * print the collection matrix
 */
void CriteriaWLTMMC::printCollectMat(const char* fileName  //!< file name
  ) {
  // determine if print reweighted or actual lnpi
  std::string colType("");
  if (printRW_) {
    colType.assign("rw");
  }

  // append state conditions on file name
  std::ostringstream colMatName;
  colMatName << fileName << colType;

  writeRestart(colMatName.str().c_str());

  // if using growth expanded ensemble with nmolstages,
  // also print physical states
  if (mType_.compare("nmolstage") == 0) {
    const int nMin = myRound(bin2m(0)), nMax = myRound(bin2m(nBin_ - 1));
    CriteriaWLTMMC cNoGrow(beta_, activ_, "nmol", nMin - 0.5, nMax + 0.5,
                           nMax - nMin + 1);
    int binNoGrow = 0;
    for (int bin = 0; bin < nBin_; ++bin) {
      if (fabs(bin2m(bin) - cNoGrow.bin2m(binNoGrow)) < 1e-12) {
        cNoGrow.lnPI_[binNoGrow] = lnPI_[bin];
        cNoGrow.pe_[binNoGrow] = pe_[bin];
        ++binNoGrow;
      }
    }
    colMatName.str("");
    colMatName << fileName << "nogrow";
    cNoGrow.printCollectMat(colMatName.str().c_str());
  }
}

/**
 * read the collection matrix and lnPI
 */
void CriteriaWLTMMC::readCollectMat(const char* fileName  //!< file name
  ) {
  std::ifstream fs(fileName);
  std::string line;
  const int nLines = numLines(fileName);

  // skip header lines
  if (nLines > nBin_) {
     for (int i = 0; i < nLines - nBin_; ++i) getline(fs, line);
  }

  for (int i = 0; i < nBin_; ++i) {
    double tmp;
    fs >> tmp >> lnPI_[i];
    getline(fs, line);
  }
  lnPInorm();
}

/**
 * zero all statistics and accumulators
 */
void CriteriaWLTMMC::zeroStat() {
  printRW_ = false;
  h_.clear();
  h_.resize(nBin_, 0);
  lnPI_.clear();
  lnPI_.resize(nBin_, -1.);
  pe_.clear();
  pe_.resize(nBin_);
  prefilColMat(0);
  // prefilColMat(std::numeric_limits<double>::min());
  lnpi2pressure_.clear();
  verbose_ = 0;
  nTunnels_ = 0;
  nTunnelPrev_ = -1;
  wlFlat_ = 0;
  nSweep_ = 0;
  lnf_ = 1.;
}

/**
 * convert collection matrix to lnPI
 */
void CriteriaWLTMMC::c2lnPI(
  const vector<vector<long double> > &col,   //!< collection matrix
  vector<long double> *lnpiPtr) {   // log of probability distribution function
  vector<long double>& lnpi = *lnpiPtr;
  double lnPIprev = 0.;
  lnpi[0] = lnPIprev;
  double p01, p10, cSum1;
  double cSum0 = std::accumulate(col[0].begin(), col[0].end(), 0.);
  for (int i = 1; i < nBin_; ++i) {
    if (cSum0 == 0) {
      lnpi[i] = lnPIprev;
      cSum0 = std::accumulate(col[i].begin(), col[i].end(), 0.);
    } else {
      if (cTripleBanded_) {
        p01 = col[i-1][2] / cSum0;
      } else {
        p01 = col[i-1][i] / cSum0;
      }
      cSum1 = std::accumulate(col[i].begin(), col[i].end(), 0.);
      cSum0 = cSum1;
      if (cSum1 == 0) {
        lnpi[i] = lnPIprev;
      } else {
        if (cTripleBanded_) {
          p10 = col[i][0] / cSum1;
        } else {
          p10 = col[i][i-1] / cSum1;
        }
        if (p10 == 0) {
          lnpi[i] = lnPIprev;
        } else {
          lnpi[i] = lnPIprev + log(p01/p10);
          lnPIprev = lnpi[i];
        }
      }
    }
  }
  lnPInorm(lnpiPtr);
}

/**
 * update macrostate probability distribution with collection matrix
 *  sum multiple collection matrices
 */
void CriteriaWLTMMC::lnPIupdate(
  const vector<std::shared_ptr<CriteriaWLTMMC> > &c) {
  if (tmmc_) {
    // sum all collection matrices into one
    vector<vector<long double> > C(static_cast<int>(
      C_.size()), vector<long double>(static_cast<int>(C_[0].size()), 0.));
    for (unsigned int i = 0; i < c.size(); ++i) {
      for (unsigned int j = 0; j < C.size(); ++j) {
        for (unsigned int k = 0; k < C[j].size(); ++k) {
          C[j][k] = c[i]->C()[j][k];
        }
      }
    }

    // update lnpi with total collection matrix
    c2lnPI(C, &lnPI_);
  }
}

/**
 * function to reweight lnPI to different value of activity
 */
void CriteriaWLTMMC::lnPIrw(
  const double activrw) {    //!< activity at reweighted condition
  activrw_ = activrw;
  // cout << "activ " << activ_ << " rw " << activrw_ << endl;
  lnPIrw_.resize(lnPI_.size());
  for (unsigned int i = 0; i < lnPI_.size(); ++i) {
    // lna = beta*mu-3lnLambda, larw - lna = beta*(murw - mu)
    const double fac = (log(activrw_) - log(activ_))*bin2m(i);
    if (fac == 0) {
      lnPIrw_[i] = lnPI_[i];
    } else {
      lnPIrw_[i] = lnPI_[i]+fac;
    }
  }
  lnPInorm(&lnPIrw_);
}

/**
 * function to find saturation by reweight lnPI to different value of activity
 *  returns squared difference of peak heights
 */
double CriteriaWLTMMC::lnPIrwsat(
  const double activrw) {    //!< activity at reweighted condition
  lnPIrw(activrw);
  vector<CriteriaWLTMMC> c = phaseSplit(lnPIrw_);
  if (c.size() != 2) {
    double returnval = 1e13;
    // encourage multiple peaks, which typically means it is flat but noisey
    if (c.size() > 2) returnval *= 1e-2;

    // bias reweights with smallest squared difference between first and last
    // lnPI value
    returnval += pow(lnPIrw_.front() - lnPIrw_.back(), 2);
    return returnval;
  } else {
    const double sqDiff = pow(log(c[0].lnPIarea(c[0].lnPI()))
                            - log(c[1].lnPIarea(c[1].lnPI())), 2.);
    return sqDiff;
  }
}

/**
 * reweight to saturation conditions by mimizing lnPIrw2sat
 */
void CriteriaWLTMMC::findSat() {
  Golden g;
  lnPIrwsatwrapper lnpirwswrap = lnPIrwsatwrap();
  g.bracket(activ_, activ_*1.1, lnpirwswrap);
  const double activ = g.minimize(lnpirwswrap);
  lnPIrw(activ);
  // cout << activrw_ << endl;
}

/**
 * initialize the macrostate bins
 */
void CriteriaWLTMMC::initBins(
  const double mMax,  //!< minimum value of macrostate
  const double mMin,  //!< maximum value of macrostate
  const int nBin) {  //!< number of bins for macrostate
  if (mType_.compare("energy") == 0) {
    cTripleBanded_ = false;
  } else if ( (mType_.compare("nmol") == 0) ||
              (mType_.compare("nmolstage") == 0) ||
              (mType_.compare("pairOrder") == 0) ) {
    cTripleBanded_ = true;
  } else if (mType_.compare("beta") == 0) {
    cTripleBanded_ = true;
    printBeta_ = 1;
  } else if ( (mType_.compare("pressure") == 0) ||
              (mType_.compare("lnpres") == 0) ) {
    cTripleBanded_ = true;
    printPressure_ = 1;
  } else {
    ASSERT(0, "unrecognized macrostate type (" << mType_ << ")");
  }

  mMin_ = mMin;
  mMax_ = mMax;
  ASSERT(mMin_ <= mMax_, "mMin(" << mMin_ << ") >= mMax(" << mMax_ << ")");
  nBin_ = nBin;
  mBin_ = ((mMax - mMin)/nBin);
  zeroStat();
}

/**
 * return new nMolMax to resize the nmol window in WLTMMC criteria
 *  using printRW to determine whether to use lnPI or lnPIrw
 *  assumes two peaks in lnPI
 *  truncates liquid peak after lnPI drops as n increases by amout liquidDrop
 *  rounds to the nearest larger multiple of nround
 */
int CriteriaWLTMMC::nMolResizeWindow(
  const double liquidDrop,  //!< targeted drop off of liquid peak
  const int round) {        //!< round up to the nearest factor
  vector<long double> *lnPItmp;
  if (printRW()) {
    lnPItmp = &lnPIrw_;
    printRW_ = false;
  } else {
    lnPItmp = &lnPI_;
  }

  // find liquid peak
  const vector<int> max = findLocalMaxima(*lnPItmp, nSmooth_);
  ASSERT(max.size() == 2, "nMolResizeWindow assumes that the number of peaks("
    << max.size() << ") = 2");
  const long double lpeak = lnPItmp->at(max.back());

  // as n increases, look for change in peak height < liquidDrop
  bool dropFound = false;
  double mx = -1;
  for (int m = max[1] + 1; m < nBin(); ++m) {
    if ( (dropFound == false) && (lnPItmp->at(m) - lpeak < liquidDrop) ) {
      dropFound = true;
      mx = bin2m(m);
    }
  }

  // expand range if it is not large enough to find the requested drop
  if (dropFound == false) {
    mx = mMax() + 2*round;
  }

  // round mx
  const int nmx = static_cast<int>(mx) - (static_cast<int>(mx) % round) + round;
  // cout << "mx " << mx << " nmx " << nmx << endl;

  return nmx;
}

/**
 * obtain pressure isotherm from lnPI
 */
void CriteriaWLTMMC::lnPIpressureIso(const double vol   //!< volume
  ) {
  lnpi2pressure_.clear();
  lnpi2pressure_.resize(nBin_);
  pressureVec_.clear();
  pressureVec_.resize(nBin_);
  volume_ = vol;
  ASSERT(fabs(bin2m(0)) < doubleTolerance, "pressure computation requires"
    << "simulation at N=0, however, bin2m(0) = " << bin2m(0) << ")");

  // first, scan for conditions where there is only one (meta)stable phase
  nMolPeakPhase_ = 0;
  double lnactivGuess = 0.1*log(activ_);
  for (int i = 1; i < nBin_-1; ++i) {
    findPeak(bin2m(i), lnactivGuess);
    lnactivGuess = log(activrw_);
    lnpi2pressure_[i] = lnPIrwpressureOnePhase(vol);
  }
}

/**
 * reweight to obtain peak maxima at given nmol
 */
void CriteriaWLTMMC::findPeak(
  const double nMolPeak,  //!< target for peak location
  const double lnactivGuess) {   //!< first guess for activity
  nMolPeak_ = nMolPeak;
  Golden g;
  lnPIrwnmxwrapper lnpirwwrap = lnPIrwnmxwrap();
  g.bracket(lnactivGuess, 1.05*lnactivGuess, lnpirwwrap);
  const double lnactiv = g.minimize(lnpirwwrap);
  lnPIrw(exp(lnactiv));
  // cout << "target " << nMolPeak << " activ " << lnactiv << endl;
}

/**
 * find peak location by reweighting and minimizing squared different in peak
 * location
 */
double CriteriaWLTMMC::lnPIrwnmx(
  const double lnactivrw) {    //!< activity at reweighted condition
  lnPIrw(exp(lnactivrw));
  return pow(lnPIaverage(lnPIrw_) - nMolPeak_, 2);
}

/**
 * split phases among vector of criteria
 */
vector<CriteriaWLTMMC> CriteriaWLTMMC::phaseSplit(
  const vector<long double> &lnPI) {
  vector<int> min = lnPIphaseBoundary(lnPI);
  const int nPhases = lnPInumPhases(lnPI);
  // cout << "nPhases " << nPhases << endl;

  // create vector of criteria with correct macrostate range
  vector<CriteriaWLTMMC> c;
  if (nPhases == 1) {
    c.push_back(*this);
  } else {
    for (int i = 0; i < nPhases; ++i) {
      double mmin, mmax;
      if (i == 0) {
        mmin = mMin_;
      } else {
        mmin = bin2m(min[i - 1]) + 0.5*mBin_;
      }
      if (i == static_cast<int>(min.size())) {
        mmax = mMax_;
      } else {
        mmax = bin2m(min[i]) + 0.5*mBin_;
      }
      // cout << "mmin " << mmin << " mmax " << mmax << endl;
      int bins;
      if (i == 0) {
        if (nPhases == 1) {
          bins = nBin_;
        } else {
          bins = min[i] + 1;
        }
      } else if (i == static_cast<int>(min.size())) {
        bins = nBin_ - min.back() - 1;
      } else {
        bins = min[i] - min[i-1];
      }
      CriteriaWLTMMC ctmp(beta_, activ_, mType_.c_str(), mmin, mmax, bins);
      c.push_back(ctmp);
    }
  }

  // put the lnPI into the vector of criteria
  // cout << "put the lnPI into a vector of criteria" << endl;
  int iPhase = 0, index = 0;
  for (unsigned int i = 0; i < lnPI.size(); ++i) {
    c[iPhase].setlnPI_(index, lnPI[i]);
    ++index;
    if (iPhase < static_cast<int>(min.size())) {
      if (static_cast<int>(i) == min[iPhase]) {
        ++iPhase;
        index = 0;
      }
    }
  }
  return c;
}

/**
 * paste together windows of lnPI to print final version
 */
void CriteriaWLTMMC::printCollectMat(
  const char* fileName,    //!< file to print aggregate lnPI
  const vector<CriteriaWLTMMC*> c) {   //!< vector of criteria with lnPI
  mout_("warning", std::ostringstream().flush()
    << "printCollectMat with vector of criteria is depreciated");

  int nWindow = static_cast<int>(c.size());
  // check if a collection matrix has not been instantiated yet
  for (int w = 0; w < nWindow; ++w) {
    if (c[w] == NULL) nWindow = -1;
  }

  if (nWindow == 0) {
    ASSERT(0, "cannot print Collection Matrix of null criteria class");
  } else if (nWindow == 1) {
    c[0]->printCollectMat(fileName);
  } else if (nWindow == -1) {
    // do nothing, not all collect mat are ready yet (beginning of run)
  } else {
    vector<double> lnPI;
    vector<double> lnPIx;
    vector<double> lnPIcut;
    vector<double> lnPIcutx;
    vector<double> pecut;
    vector<double> pecutstd;
    for (unsigned int i = 0; i < c.front()->lnPI().size(); ++i) {
      lnPI.push_back(c.front()->lnPI()[i]);
      lnPIx.push_back(c.front()->bin2m(i));
      lnPIcut.push_back(c.front()->lnPI()[i]);
      lnPIcutx.push_back(c.front()->bin2m(i));
      pecut.push_back(c.front()->pe()[i].average());
      pecutstd.push_back(c.front()->pe()[i].stdev());
    }
    // for each window!=0, find where first macrostate overlaps previous window
    // and shift entire window by that difference
    double shift = 0;
    for (int w = 1; w < nWindow; ++w) {
      // found overlap shift constant
      bool foundOverlap = false;
      int nOverlap = 0;
      int overlapingBins = 0;
      for (unsigned int i = 0; i < c[w-1]->lnPI().size(); ++i) {
        if (fabs(c[w-1]->bin2m(i) - c[w]->bin2m(0)) < 1e-12) {
          if (foundOverlap == false) {
            foundOverlap = true;
            shift += c[w-1]->lnPI()[i] - c[w]->lnPI()[0];
            nOverlap = c[w-1]->nBin() - i;
            overlapingBins = c[w-1]->nBin() - i;
          }
        }
      }

      const int npop = static_cast<int>((overlapingBins-1)/2);
      nOverlap -= npop;
      for (int i = 0; i < npop; ++i) {
        lnPIcut.pop_back();
        lnPIcutx.pop_back();
        pecut.pop_back();
        pecutstd.pop_back();
      }

      if (foundOverlap) {
        for (unsigned int i = 0; i < c[w]->lnPI().size(); ++i) {
          lnPI.push_back(c[w]->lnPI()[i] + shift);
          lnPIx.push_back(c[w]->bin2m(i));
          if (static_cast<int>(i) >= nOverlap) {
            lnPIcut.push_back(c[w]->lnPI()[i] + shift);
            lnPIcutx.push_back(c[w]->bin2m(i));
            pecut.push_back(c[w]->pe()[i].average());
            pecutstd.push_back(c[w]->pe()[i].stdev());
          }
        }
      }
    }

    // normalize lnPIs
    lnPInorm(&lnPI);
    lnPInorm(&lnPIcut);

    // print
    // append state conditions on file name
    std::ostringstream colMatName;
    colMatName << fileName;

    fileBackUp(colMatName.str().c_str());
    {
      std::ofstream file(colMatName.str().c_str(), std::ios::trunc);
      for (unsigned int i = 0; i < lnPI.size(); ++i) {
        file << lnPIx[i] << " " << lnPI[i] << endl;
      }
    }
    fileBackUp(colMatName.str().c_str());
    lnPIaggre_.clear();
    {
      std::ofstream file(colMatName.str().c_str(), std::ios::trunc);
      for (unsigned int i = 0; i < lnPIcut.size(); ++i) {
        file << lnPIcutx[i] << " " << lnPIcut[i] << " " << pecut[i] << " "
             << pecutstd[i] << endl;
        lnPIaggre_.push_back(lnPIcut[i]);
      }
    }
  }
}

/**
 * replace lnPI
 */
void CriteriaWLTMMC::lnPIreplace(const vector<long double> &lnPI) {
  // update lnPI of c_ to aggregate of windows
  ASSERT(lnPI_.size() == lnPI.size(), "lnPI(size=" << lnPI.size()
    << ") doesn't matc lnPI(size=" << lnPI_.size() << ")");
  for (unsigned int i = 0; i < lnPI_.size(); ++i) {
    lnPI_[i] = lnPI[i];
  }
}

/**
 * obtain energy isotherm from lnPI
 */
void CriteriaWLTMMC::lnPIenergyIso() {
  vector<double> pe;
  for (int i = 0; i < nBin_; ++i) {
    pe.push_back(pe_[i].average());
  }
  peMUVT_ = lnPIgc2can(pe);
}

/**
 * read lnPI and energy
 */
void CriteriaWLTMMC::readlnPIEner(const char* fileName  //!< file name
  ) {
  std::ifstream fs(fileName);
  std::string line;
  const int nLines = numLines(fileName);

  // skip header lines
  if (nLines > nBin_) {
     for (int i = 0; i < nLines - nBin_; ++i) getline(fs, line);
  }

  // // skip all lines beginning with the character "#"
  // skipCharsInFile('#', fs);

  for (int i = 0; i < nBin_; ++i) {
    double tmp;
    fs >> tmp >> lnPI_[i] >> tmp;
    pe_[i].accumulate(tmp);
    getline(fs, line);
  }
  lnPInorm();
}

/**
 * splice windows together
 */
void CriteriaWLTMMC::spliceWindows(const vector<CriteriaWLTMMC*> c) {
  int nWindow = static_cast<int>(c.size());
  // check if a collection matrix has not been instantiated yet
  for (int w = 0; w < nWindow; ++w) {
    if (c[w] == NULL) nWindow = -1;
  }

  ASSERT(nWindow != 0, "cannot print Collection Matrix of null criteria class");
  if (nWindow == 1) {
    mout_("warning", std::ostringstream().flush()
      << "no reason to splice only one window");
  } else if (nWindow == -1) {
    // do nothing, not all collect mat are ready yet (beginning of run)
  } else {
    nSweep_ = minNSweep(c);
    int bin = 0;
    for (unsigned int i = 0; i < c.front()->lnPI().size(); ++i) {
      lnPI_[bin] = c.front()->lnPI()[i];
      pe_[bin] = c.front()->pe()[i];
      ++bin;
    }
    // for each window!=0, find where first macrostate overlaps previous window
    // and shift entire window by that difference
    double shift = 0;
    for (int w = 1; w < nWindow; ++w) {
      // found overlap shift constant
      bool foundOverlap = false;
      int nOverlap = 0;
      int overlapingBins = 0;
      for (unsigned int i = 0; i < c[w-1]->lnPI().size(); ++i) {
        if (fabs(c[w-1]->bin2m(i) - c[w]->bin2m(0)) < 1e-12) {
          if (foundOverlap == false) {
            foundOverlap = true;
            shift += c[w-1]->lnPI()[i] - c[w]->lnPI()[0];
            nOverlap = c[w-1]->nBin() - i;
            overlapingBins = c[w-1]->nBin() - i;
          }
        }
      }

      const int npop = static_cast<int>((overlapingBins-1)/2);
      nOverlap -= npop;
      bin -= npop;

      if (foundOverlap) {
        for (unsigned int i = 0; i < c[w]->lnPI().size(); ++i) {
          if (static_cast<int>(i) >= nOverlap) {
            lnPI_[bin] = c[w]->lnPI()[i] + shift;
            pe_[bin] = c[w]->pe()[i];
            ++bin;
          }
        }
      }
    }

    // normalize lnPI
    lnPInorm();
  }
}

/**
 * minimum number of sweeps in all windows
 */
int CriteriaWLTMMC::minNSweep(const vector<CriteriaWLTMMC*> c) {
  if (c.size() == 0) {
    return 0;
  } else {
    int minSweep = c.front()->nSweep();

    // check if a collection matrix has not been instantiated yet
    int nWindow = static_cast<int>(c.size());
    for (int w = 0; w < nWindow; ++w) {
      if (c[w] == NULL) nWindow = -1;
    }

    if (nWindow > 0) {
      for (unsigned int i = 1; i < c.size(); ++i) {
        if (c[i]->nSweep() < minSweep) minSweep = c[i]->nSweep();
      }
      return minSweep;
    } else {
      return 0;
    }
  }
}

/**
 * minimum number of sweeps in all windows
 */
int CriteriaWLTMMC::minNwlFlat(const vector<CriteriaWLTMMC*> c) {
  if (c.size() == 0) {
    return 0;
  } else {
    int minWLFlat = c.front()->wlFlat();

    // check if a collection matrix has not been instantiated yet
    int nWindow = static_cast<int>(c.size());
    for (int w = 0; w < nWindow; ++w) {
      if (c[w] == NULL) nWindow = -1;
    }

    if (nWindow > 0) {
      for (unsigned int i = 1; i < c.size(); ++i) {
        if (c[i]->wlFlat() < minWLFlat) minWLFlat = c[i]->wlFlat();
      }
      return minWLFlat;
    } else {
      return 0;
    }
  }
}

/**
 * prefill collection matrix with a constant
 */
void CriteriaWLTMMC::prefilColMat(const long double constant) {
  C_.clear();
  if (cTripleBanded_) {
    C_.resize(nBin_, vector<long double>(3, constant));
  } else {
    C_.resize(nBin_, vector<long double>(nBin_, constant));
  }
}

/**
 * write restart file
 */
void CriteriaWLTMMC::writeRestart(const char* fileName) {
  // determine if print reweighted or actual lnpi
  vector<long double> lnPItmp;
  vector<long double> lnPIwlcomp;
  if (printRW_) {
    lnPItmp = lnPIrw_;
  } else {
    lnPItmp = lnPI_;
  }
  lnPInorm(&lnPItmp);

  if (collect_ && !tmmc_) {
    lnPIwlcomp.resize(lnPItmp.size());
    c2lnPI(C_, &lnPIwlcomp);
  }

  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);

  std::streamsize ss = cout.precision();

  // state conditions and run parameters
  file << "# collect " << collect_ << endl
       << "# tmmc " << tmmc_ << endl
       << "# nSweep " << nSweep_ << endl
       << "# wlFlat " << wlFlat_ << endl
       << "# lnf " << lnf_ << endl
       << "# gwlmod " << g_ << endl
       << "# mType " << mType_ << endl
       << std::setprecision(std::numeric_limits<long double>::digits10+2)
       << "# mMin " << mMin_ << endl
       << "# mMax " << mMax_ << endl
       << std::setprecision(ss)
       << "# nBin " << nBin_ << endl
       << "# lnfCollect " << lnfCollect_ << endl
       << "# lnfTMMC " << lnfTMMC_ << endl
       << "# wlFlatFactor " << wlFlatFactor_ << endl
       << "# nSweepVisPerBin " << nSweepVisPerBin_ << endl;

  // header
  file << "# macrostate(" << mType_ << ") lnPi(m) ";
  if (lnpi2pressure_.size() == C_.size()) file << "rho pressure ";
  file << "pe pe_stdev ";
  if (peMUVT_.size() == C_.size()) file << "peMUVT ";
  if (cTripleBanded_) file << "colMat(m-1) colMat(m) colMat(m+1) ";
  if (collect_ && !tmmc_) file << "lnPIwlcomp ";
  file << "h peNvalues peSum peSumSq" << endl;

  for (unsigned int i = 0; i < C_.size(); ++i) {
    file << bin2m(i) << " ";
    file << std::setprecision(std::numeric_limits<long double>::digits10+2)
         << lnPItmp[i] << " ";
    file << std::setprecision(ss);
    if (lnpi2pressure_.size() == C_.size()) {
      file << bin2m(i)/volume_ << " " << lnpi2pressure_[i] << " ";
    }
    file << pe_[i].average() << " " << pe_[i].stdev() << " ";
    if (peMUVT_.size() == C_.size()) {
      file << peMUVT_[i] << " ";
    }
    for (unsigned int j = 0; j < C_[i].size(); ++j) {
      file << std::setprecision(std::numeric_limits<long double>::digits10+2)
           << C_[i][j] << " ";
    }
    file << std::setprecision(ss);
    if (collect_ && !tmmc_) file << lnPIwlcomp[i] << " ";
    file << h_[i] << " ";
    file << std::setprecision(std::numeric_limits<long double>::digits10+2)
         << pe_[i].nValues() << " "
         << pe_[i].sum() << " "
         << pe_[i].sumSq() << " "
         << pe_[i].blockStdev() << endl;
    file << std::setprecision(ss);
  }
  printRW_ = false;
}

/**
 * reweight to find difference in peaks
 */
void CriteriaWLTMMC::peakDiff(const int iPeak, const int jPeak) {
  if (iPeak == jPeak) {}  // remove warning for unused parameters
  nMolPeakPhase_ = 0;
  double lnactivGuess = 0.1*log(activ_);
  for (int ibin = 1; ibin < nBin_-1; ++ibin) {
    findPeak(bin2m(ibin), lnactivGuess);
    lnactivGuess = log(activrw_);
    vector<int> max = findLocalMaxima(lnPIrw_, 3);
    max.resize(10);
    cout << ibin << " ";
    for (unsigned int i = 0; i < max.size(); ++i) {
      cout << lnPIrw_[max[i]] - lnPIrw_[max[0]] << " ";
    }
    cout << endl;
  }
}

/**
 * read lnPI, energy and collection matrix
 */
void CriteriaWLTMMC::readlnPIEnerCol(const char* fileName) {  //!< file name
  std::ifstream fs(fileName);
  std::string line;
  const int nLines = numLines(fileName);

  // skip header lines
  if (nLines > nBin_) {
     for (int i = 0; i < nLines - nBin_; ++i) getline(fs, line);
  }

  for (int i = 0; i < nBin_; ++i) {
    vector<double> tmp(6);
    fs >> tmp[0] >> lnPI_[i] >> tmp[1] >> tmp[2] >> tmp[3] >> tmp[4] >> tmp[5];
    pe_[i].accumulate(tmp[1]);
    if (cTripleBanded_) {
      C_[i][0] = tmp[3];
      C_[i][1] = tmp[4];
      C_[i][2] = tmp[5];
    } else {
      ASSERT(0,
        "reading non-triple banded collection matrix is not implemented");
    }
    if (collect_ && !tmmc_) {
      fs >> tmp[0];
    }
    int h;
    long long nValues;
    long double sum, sumSq;
    fs >> h >> nValues >> sum >> sumSq >> tmp[0];
    h_[i] = h;
    Accumulator pe(nValues, sum, sumSq);
    pe_[i] = pe;
    getline(fs, line);
  }
  lnPInorm();
}

/**
 * obtain grand canonical ensemble average from canonical ensemble average
 */
vector<double> CriteriaWLTMMC::lnPIgc2can(const vector<double> data) {
  vector<double> returnData(data.size());
  nMolPeakPhase_ = 0;
  double lnactivGuess = 0.1*log(activ_);
  for (int i = 1; i < nBin_-1; ++i) {
    findPeak(bin2m(i), lnactivGuess);
    lnactivGuess = log(activrw_);
    returnData[i] = lnPIrwdata(data);
  }
  return returnData;
}

/**
 * obtain grand canonical ensemble average from canonical ensemble average
 *  input file has macrostate as first column, data as second column
 */
void CriteriaWLTMMC::lnPIgc2can(const char* fileNameIn,
  const char* fileNameOut) {
  // read canonical ensemble data from file
  vector<double> dataNVT(nBin_);
  std::ifstream fileIn(fileNameIn);

  // skip all lines beginning with the character "#"
  // skipCharsInFile('#', fileIn);

  string line;
  int index;
  double data;
  while (!fileIn.eof()) {
    fileIn >> index >> data;
    cout << "index " << index << " data " << data << endl;
    if ( (index < 0) || (index > nBin_ - 1) ) {
      mout_("warning", std::ostringstream().flush()
        << "data from file(" << fileNameIn << ") is out of bounds, index("
        << index << ") > nBin(" << nBin_ << ")");
      getline(fileIn, line);
    } else {
      dataNVT[index] = data;
    }
    getline(fileIn, line);
  }

  // reweight to grand canonical ensemble
  const vector<double> dataMUVT = lnPIgc2can(dataNVT);

  // output grand canonical ensemble data
  std::ofstream fileOut(fileNameOut);
  for (unsigned int i = 0; i < dataMUVT.size(); ++i) {
    fileOut << i << " " << dataMUVT[i] << endl;
  }

  vector<int> max = findLocalMaxima(dataMUVT, 1);
  cout << dataMUVT[max.front()] << endl;
}

/**
 * return heat capacity for each bin
 */
vector<double> CriteriaWLTMMC::heatCapacity() {
  vector<double> cv;
  for (int bin = 0; bin < nBin(); ++bin) {
    const double nValues = pe_[bin].nValues(),
      sum   = pe_[bin].sum(),
      sumSq = pe_[bin].sumSq();
    cv.push_back(((sumSq/nValues)-pow(sum/nValues, 2.))*pow(beta_, 2.) );
  }
  return cv;
}

/**
 * return energy fluctuations for each bin
 */
vector<double> CriteriaWLTMMC::fluct() {
  vector<double> cv;
  for (int bin = 0; bin < nBin(); ++bin) {
    const double nValues = pe_[bin].nValues(),
      sum   = pe_[bin].sum(),
      sumSq = pe_[bin].sumSq();
    cv.push_back(((sumSq/nValues)-pow(sum/nValues, 2.)));
  }
  return cv;
}

/**
 * initialize TMMC
 */
void CriteriaWLTMMC::tmmcInit() {
  ASSERT(collect_, "must also fill collection matrix if running tmmc");
  tmmc_ = true;
  std::fill(h_.begin(), h_.end(), 0);
}

}  // namespace feasst

