/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef CRITERIA_WLTMMC_H_
#define CRITERIA_WLTMMC_H_

#include <math.h>
#include "./functions.h"
#include "./criteria.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Wang-Landau (WL) and Transition Matrix (TM) Monte Carlo acceptance criteria.
 * http://dx.doi.org/10.1063/1.1572463
 * http://dx.doi.org/10.1063/1.4884124
 */
class CriteriaWLTMMC : public Criteria {
 public:
  /// Constructor.
  CriteriaWLTMMC(const double beta,
    /**
     * allowed string key pairs (e.g., dictionary)
     *
     * mType : type of macrostate
     *  - energy : potential energy of system
     *  - nmol : number of molecules
     *  - nmol0 : number of molecules of the first type
     *  - nmolstage : number of molecules with growth expanded ensemble
     *  - pairOrder : order parameter defined by the Pair class
     *  - beta : temperature expanded ensemble
     *  - pressure : thermodynamic pressure
     *  - lnpres : logarithmic pressure
     *
     * mMax : maximum floating point value of macrostate.
     * mMin : minimum floating point value of macrostate.
     *
     * mMaxCenter : maximum center bin floating point value of macrostate,
     *   if mMax not provided.
     * mMinCenter : minimum center bin floating point value of macrostate.
     *
     * nBin : number of equal sized bins for the 1D macrostate range.
     *
     * nMax : maximum integer value of macrostate, if mMax(Center) not provided.
     *   Warning: do not provide nBin and nMax together as arguments.
     * nMin : minimum integer value of macrostate (default: 0).
     */
    const argtype &args);

  // Constructor. Arguments are as described below.
  // HWH Depreciate in favor of above
  CriteriaWLTMMC(const double beta, const double activ, const char* mType,
                 const double mMin, const double mMax, const int nBin);

  /** Return macrostate type or name. Value types are as follows:
   *  "energy" for potential energy,
   *  "nmol" for number of molecules,
   *  "nmol0" for number of molecules of the first initialized type,
   *  "nmolstage", for number of molecules with growth expanded ensemble,
   *  "pairOrder" for order parameter defined by the Pair class,
   *  "beta" for inverse temperature,
   *  "pressure" for the thermodynamic pressure, and
   *  "lnpres" for logarithmic pressure. */
  string mType() const { return mType_; }

  /// Return minimum value of macrostate.
  double mMin() const { return mMin_; }

  /// Return maximum value of macrostate.
  double mMax() const { return mMax_; }

  /// Return number of macrostate bins between mMin and mMax.
  int nBin() const { return nBin_; }

  /// Return size of macrostate bins, which are constant.
  double mBin() const { return mBin_; }

  /// Alternative constructor. In this case, integer bin size is assumed.
  CriteriaWLTMMC(const double beta, const double activ, const char* mType,
                 const int nMin, const int nMax);

  /// Initialize the macrostate bins, which are constant in size.
  void initBins(const double mMax, const double mMin, const int nBin);

  /// Return whether to accept (1) or reject (0) the proposed trial.
  int accept(const double pMet, const double peNew, const char* moveType,
             const int reject);

  /// Store macrostate variables of old configuration.
  void store(Pair* pair);

  /// Return value at center of bin, given bin.
  double bin2m(const int bin) const { return mMin_ + (bin + 0.5)*mBin_; }

  /// Return bin, given value.
  int bin(const double m) const { return feasstRound((m - bin2m(0))/mBin_); }

  // Functions for Wang-Landau

  /// Set the Wang-Landau modification factor. Otherwise, default value.
  void setg(const double g = 0.5) { g_ = g; }

  /** Set the Wang-Landau flatness threshold factor for flatness checks:
   *  min(histogram) >= factor*average(histogram) */
  void setWLFlatFactor(const double factor = 0.8) { wlFlatFactor_ = factor; }

  /** Return 1 if flatness criteria is met. Reset histogram and reduce update
   *  factor(lnf) by Wang-Landau modification factor. */
  void flatCheck(const int force = 0  //!< force flatness true if 1
    );

  /// Return the Wang-Landau update factor.
  double lnf() const { return lnf_; }

  /// Return the number of Wang-Landau flatness criteria that have been met.
  int wlFlat() const { return wlFlat_; }

  // Functions for transition-matrix

  /** Initialize the collection matrix, but not necessarily the use of the
   *  resulting macrostate probability for selecting trials (e.g., TMMC). */
  void collectInit() {collect_ = true; }

  /** Initialize updating of the collection matrix for when the value of the
   *  Wang-Landau update factor reaches lnfCollect. */
  void collectInit(const double lnfCollect) {lnfCollect_ = lnfCollect; }

  /** Initialize updating of the collection matrix for when the value of the
   *  number of Wang-Landau flatness checks reaches wlFlat. */
  void collectInit(const int wlFlat) {
    lnfCollect_ = (pow(g_, wlFlat) + pow(g_, wlFlat - 1))/2.0;
  }

  /// Initialize Transition-Matrix Monte Carlo.
  void tmmcInit();

  /** Initialize Transition-Matrix Monte Carlo for when the value of the
   *  Wang-Landau update factor reaches lnfTMMC. */
  void tmmcInit(const double lnfTMMC) {lnfTMMC_ = lnfTMMC; }

  /** Initialize Transition-Matrix Monte Carlo for when the value of the
   *  number of Wang-Landau flatness checks reaches wlFlat. */
  void tmmcInit(const int wlFlat) {
    lnfTMMC_ = (pow(g_, wlFlat) + pow(g_, wlFlat - 1))/2.0;
  }

  /// Update macrostate probability distribution with collection matrix.
  void lnPIupdate();

  /** Update macrostate probability distribution, lnPI with multiple collection
   *  matrices. */
  void lnPIupdate(const vector<std::shared_ptr<CriteriaWLTMMC> > &c);

  /// Convert collection matrix, col, to probability distribution, lnPI.
  void c2lnPI(const vector<vector<long double> > &col,
              vector<long double> *lnpiPtr);

  /// Return area for a given lnPI.
  template<class T>
  double lnPIarea(const vector<T> &lnPI) const {
    double sum = 0;
    for (unsigned int i = 0; i < lnPI.size(); ++i) sum += exp(lnPI[i]);
    return sum;
  }

  /// Return area of the current lnPI.
  double lnPIarea() const { return lnPIarea(lnPI_); }

  /// Return area of the reweighted lnPI.
  double lnPIrwarea() const { return lnPIarea(lnPIrw_); }

  /// Normalize lnPI such that Sum(Pi) = 1.
  template<class T>
  void lnPInorm(vector<T> *lnPI) {
    //  to avoid exp overflow, shift by const highest value
    if (lnPI->size() > 0) {
      const T shift = *std::max_element(lnPI->begin(), lnPI->end());
      for (unsigned int i = 0; i < lnPI->size(); ++i) (*lnPI)[i] -= shift;
      const T lns = log(lnPIarea(*lnPI));
      for (unsigned int i = 0; i < lnPI->size(); ++i) (*lnPI)[i] -= lns;
    }
  }

  /// Normalize the current lnPI.
  void lnPInorm() { lnPInorm(&lnPI_); }

  /// Return the ensemble averaged macrostate from lnPI.
  template<class T>
  double lnPIaverage(const vector<T> &lnPI) {
    double av = 0;
    for (unsigned int m = 0; m < lnPI.size(); ++m) {
      av += bin2m(m) * exp(lnPI[m]);
    }
    return av/lnPIarea(lnPI);
  }

  /// Return the ensemble average macrostate of the current lnPI.
  double lnPIaverage() { return lnPIaverage(lnPI_); }

  /// Return the ensemble average macrostate of the reweighted lnPI.
  double lnPIrwaverage() { return lnPIaverage(lnPIrw_); }

  /** Set the index of bin of the macrostate that is the boundary between
   *  phasesset phase. */
  void setPhaseBoundary(const int index) { phaseBoundary_ = index; }

  /* Return phase boundaries, e.g. lnPI minima that are not first and last
   * points. */
  template<class T>
  vector<int> lnPIphaseBoundary(const vector<T> &lnPI) {
    vector<int> min;
    if (phaseBoundary_ != 0) {
      min.push_back(phaseBoundary_);
    } else {
      min = findLocalMinima(lnPI, nSmooth_);
      if (min.size() != 0) {
        if (min.front() == 0) min.erase(min.begin());
        if (min.size() != 0) {
          if (min.back() == nBin_ - 1) min.pop_back();
        }
        const int minsize = static_cast<int>(min.size());
        if (minsize > 1) {
          min[0] = min.back();
          for (int i = 1; i < minsize; ++i) {
            min.pop_back();
          }
        }
      }
    }
    return min;
  }

  // Return phase boundaries for current macrostate.
  vector<int> lnPIphaseBoundary() { return lnPIphaseBoundary(lnPI_); }

  // Return phase boundaries for reweighted macrostate.
  vector<int> lnPIrwphaseBoundary() { return lnPIphaseBoundary(lnPIrw_); }

  // Return the number of phases.
  template<class T>
  int lnPInumPhases(const vector<T> &lnPI) {
    vector<int> min = lnPIphaseBoundary(lnPI);
    return 1 + static_cast<int>(min.size());
  }

  // Return the number of phases of the current lnPI.
  int lnPInumPhases() { return lnPInumPhases(lnPI_); }

  // Return the number of phases of the reweighted lnPI.
  int lnPIrwnumPhases() { return lnPInumPhases(lnPIrw_); }

  // Split phases among vector of criteria
  vector<CriteriaWLTMMC> phaseSplit(const vector<long double> &lnPI);

  // split lnPI among phases
  template<class T>
  vector<vector<T> > lnPIsplit(const vector<T> &lnPI) {
    vector<int> min = lnPIphaseBoundary(lnPI);
    const int nPhases = lnPInumPhases(lnPI);
    vector<vector<T> > lnPIs(nPhases);
    int iPhase = 0;
    for (unsigned int i = 0; i < lnPI.size(); ++i) {
      if (iPhase < static_cast<int>(min.size())) {
        if (i == min[iPhase]) {
          ++iPhase;
        }
      }
      lnPIs[iPhase].push_back(lnPI[i]);
    }
    return lnPIs;
  }


  // Print collection matrix and lnPI
  void printCollectMat(const char* fileName);

  // paste together windows of collection matrices to print final version
  void printCollectMat(const char* fileName, const vector<CriteriaWLTMMC*> c);
  void printCollectMat(const char* fileName,
    const vector<shared_ptr<CriteriaWLTMMC> > c) {
    return printCollectMat(fileName, shrPtr2Raw(c));
  }

  // read collection matrix
  void readCollectMat(const char* fileName);
  void readlnPIEner(const char* fileName);
  void readlnPIEnerCol(const char* fileName);

  /// Zero all statistics and accumulators.
  void zeroStat();

  /// Reweight lnPI to different value of activity.
  void lnPIrw(const double activrw);

  /** Reweight to saturation conditions my minizing the differences in peak
   *  heights. */
  void findSat();

  /** For one time, print reweighted macrostate instead of the current
   *  macrostate. */
  void printRWinit() { printRW_ = true; }

  /// Reweight to obtain peak maxima at given nmol.
  void findPeak(const double nMolPeak,  //!< target peak location
    const double lnactivGuess   //!< first guess for activity
    );

//  /// reweight to find difference in peaks
//  void peakDiff(const int iPeak, const int jPeak);

  /** return new nMolMax to resize the nmol window in WLTMMC criteria
   *  using printRW to determine whether to use lnPI or lnPIrw
   *  assumes two peaks in lnPI
   *  truncates liquid peak after lnPI drops as n increases by amout liquidDrop
   *  rounds to the nearest larger multiple of nround */
  int nMolResizeWindow(const double liquidDrop, const int round);

  /// Return ressure from lnPI.
  template<class T>
  vector<double> lnPIpressureVec(const vector<T> &lnPI, const double vol,
    const vector<CriteriaWLTMMC> &cvec) {
    vector<double> p;
    for (unsigned int i = 0; i < cvec.size(); ++i) {
      p.push_back((-lnPI.front() + log(cvec[i].lnPIarea()))/vol/beta_);
    }
    return p;
  }
  template<class T>
  double lnPIpressure(const vector<T> &lnPI, const double vol,
    const vector<CriteriaWLTMMC> &cvec) {
    vector<double> p = lnPIpressureVec(lnPI, vol, cvec);
    double maxVal = 0;
    int maxIndex = -1;
    for (unsigned int i = 0; i < cvec.size(); ++i) {
      if (cvec[i].lnPIarea() > maxVal) {
        maxIndex = i;
        maxVal = cvec[i].lnPIarea();
      }
    }
    return p[maxIndex];
  }
  double lnPIpressure(const double vol) {
    vector<CriteriaWLTMMC> cvec = phaseSplit(lnPI_);
    return lnPIpressure(lnPI_, vol, cvec);
  }
  double lnPIrwpressure(const double vol) {
    vector<CriteriaWLTMMC> cvec = phaseSplit(lnPIrw_);
    return lnPIpressure(lnPIrw_, vol, cvec);
  }
  vector<double> lnPIrwpressureVec(const double vol) {
    vector<CriteriaWLTMMC> cvec = phaseSplit(lnPIrw_);
    return lnPIpressureVec(lnPIrw_, vol, cvec);
  }
  double lnPIrwpressureOnePhase(const double vol) {
    return (-lnPIrw_.front() + log(lnPIrwarea()))/vol/beta_;
  }

  /// Obtain pressure isotherm from lnPI.
  void lnPIpressureIso(const double volume);

  /// Obtain grand canonical ensemble average from canonical ensemble average.
  vector<double> lnPIgc2can(vector<double> data);

  /**
   * obtain grand canonical ensemble average from canonical ensemble average
   *  input file has macrostate as first column, data as second column
   */
  void lnPIgc2can(const char* fileNameIn, const char* fileNameOut);

  // data from lnPI
  template<class T>
  double lnPIdata(const vector<double> data,
    const vector<T> &lnPI) {
    ASSERT(data.size() == lnPI.size(), "lnPIdata size mismatch");
    double av = 0;
    for (unsigned int m = 0; m < lnPI.size(); ++m) {
      av += data[m] * exp(lnPI[m]);
    }
    return av/lnPIarea(lnPI);
  }
  double lnPIdata(const vector<double> data) { return lnPIdata(data, lnPI_); }
  double lnPIrwdata(const vector<double> data) {
    return lnPIdata(data, lnPIrw_);
  }

  /// Obtain energy isothermo from lnPI.
  void lnPIenergyIso();

  /// Replace lnPI.
  void lnPIreplace(const vector<long double> &lnPI);

  /// Splice windows together.
  void spliceWindows(const vector<CriteriaWLTMMC*> c);

  /// Splice windows together.
  void spliceWindows(const vector<shared_ptr<CriteriaWLTMMC> > c)
    { return spliceWindows(shrPtr2Raw(c)); }

  /// Return minimum number of sweeps in all windows.
  int minNSweep(const vector<CriteriaWLTMMC*> c);

  /// Return minimum number of sweeps in all windows.
  int minNSweep(const vector<shared_ptr<CriteriaWLTMMC> > c)
    { return minNSweep(shrPtr2Raw(c)); }

  /// Return minimum number of wlFLat in all windows.
  int minNwlFlat(const vector<CriteriaWLTMMC*> c);

  /// Return minimum number of wlFLat in all windows.
  int minNwlFlat(const vector<shared_ptr<CriteriaWLTMMC> > c)
    { return minNwlFlat(shrPtr2Raw(c)); }

  /// Prefill collection matrix with a constant.
  void prefilColMat(const long double constant);

  /// Return heat capacity for each bin.
  vector<double> heatCapacity();

  /// Return energy fluctations for each bin.
  vector<double> fluct();

  /// Return current macrostate updated every call to accept() or every trial.
  int iMacro() const { return bin(mNew_); }

  /// Set the highest order of moments recorded for potential energy.
  void initMoments(const int nMoments);

  /// read-only functions for protected variables
  vector<long double> lnPI() { return lnPI_; }
  vector<long double> lnPIrw() const { return lnPIrw_; }
  vector<long double> lnPIaggre() const { return lnPIaggre_; }
  vector<vector<long double> > C() const { return C_; }
  double lastbin2m() const { return bin2m(nBin_ - 1); }
  double mOld() const { return mOld_; }
  double mNew() const { return mNew_; }
  double nTunnels() const { return nTunnels_; }
  double nSweep() const { return nSweep_; }
  double activrw() const { return activrw_; }
  bool printRW() const { return printRW_; }
  bool collect() const { return collect_; }
  bool tmmc() const { return tmmc_; }
  vector<double> lnpi2pressure() const { return lnpi2pressure_; }
  vector<Accumulator> pe() const { return pe_; }
  Accumulator pe(const int bin) { return pe_[bin]; }
  double lnfCollect() const { return lnfCollect_; }
  vector<double> peMUVT() { lnPIenergyIso(); return peMUVT_; }
  vector<shared_ptr<CriteriaWLTMMC> > crits() { return crits_; }

  /// Return estimated standard deviation of lnPI from replicas
  double lnPIstd(const double iMacro);

  /// Construct by checkpoint file.
  explicit CriteriaWLTMMC(const char* fileName);

  /// Write restart/checkpoint file.
  void writeRestart(const char* fileName);

  ~CriteriaWLTMMC() {}
  CriteriaWLTMMC* clone() const;
  shared_ptr<CriteriaWLTMMC> cloneShrPtr() const;

 protected:
  string mType_;      //!< definition of macrostate
  double mMin_;       //!< minimum value of macrostate
  double mMax_;      //!< maximum value of macrostate
  int nBin_;          //!< number of bins for macrostate
  double mBin_;       //!< macrostate bin size
  double lnf_;        //!< update factor for macrostate probability distribution
  double g_;          //!< parameter to modify update factor (0, 1) default: 1/2
  vector<int> h_;     //!< visited states histogram
  vector<long double> lnPI_;      //!< macrostate probability distribution
  double mOld_;               //!< macrostate of old configuration;
  double mNew_;               //!< macrostate of new configuration;
  double peOld_;              //!< potential energy of old configuration;
  vector<Accumulator> pe_;    //!< nvt potential energy
  vector<double> peMUVT_;     //!< muvt potential energy, U^MUVT=sum(PI(N)*U(N))
  bool cTripleBanded_;        //!< is collection matrix triple banded
  vector<vector<long double> > C_;     //!< collection matrix

  /// number of times the system moves back and forth between mMin to mMax
  int nTunnels_;

  /// previous tunnel is mMin if 0, mMax if 1, and -1 if no previous tunnels
  int nTunnelPrev_;

  /// number of times a bin must be visited (from a different bin) to count
  //  toward a sweep
  int nSweepVisPerBin_;

  /// number of times the system sweeps all bins between mMin to mMax
  int nSweep_;
  int wlFlat_;        //!< number of times Wang-Landau flatness criteria is met
  /// flatness factor for flatness criteria, min(h)>=0.8*average(h)
  double wlFlatFactor_;

  /// when lnf_ is lower than this value, begin filling of collection matrix
  double lnfCollect_;
  double lnfTMMC_;     //!< when lnf_ is lower than this value, begin using tmmc
  bool collect_;       //!< populates collection matrix, if true

  /// uses collection matrix to obtain lnPI_ (default: false)
  bool tmmc_;
  vector<long double> lnPIaggre_;  //!< aggregated lnPI from windows

  // reweighting variables
  /// reweighted macrostate probability distribution
  vector<long double> lnPIrw_;
  double activrw_;      //!< reweighted activity
  bool printRW_;      //!< print reweighted distribution instead of current lnPI
  double nMolPeak_;        //!< target for peak location
  int nMolPeakPhase_;   //!< target phase for peak location
  int nSmooth_;       //!< number of macrostates to smooth by to find min/max

  // isotherm from tmmc
  vector<double> lnpi2pressure_;       //!< pressure obtained from lnPI
  vector<vector<double> > pressureVec_;       //!< pressure obtained from lnPI
  double volume_;                 //!< volume of space, input from pressure

  void setlnPI_(const int index, const long double &val) { lnPI_[index] = val; }
  int phaseBoundary_;           //!< bin index for setting phase boundary

  /// defaults in constructor
  void defaultConstruction_();

  /// Update histograms and collection matrix for the visited macrostate.
  void mUpdate_(const int mOldBin, const int mNewBin, const double pmet,
               const int acceptFlag);

  /// Separate collection matrices for statistical error analysis
  vector<shared_ptr<CriteriaWLTMMC> > crits_;
  int separate_ = 0;

  /** Return the squared difference of peak heights after reweighting to new
   *  activity. */
  double lnPIrwsat_(const double activrw);

  /** Return the squared difference of reweighted average macrostate and
   *  the target macrostate peak, nMolPeak_. */
  double lnPIrwnmx_(const double lnactivrw);

  /// functor wrappers to pass to numerical recipe minimization algorithms
  struct lnPIrwsat_wrapper_ {
    explicit lnPIrwsat_wrapper_(CriteriaWLTMMC* this_)
    : this_(this_) {}
    double operator ( )(const double & value) {
      return this_->lnPIrwsat_(value);
    }
    CriteriaWLTMMC* this_;
  };

  lnPIrwsat_wrapper_ lnPIrwsat_wrap_() {
    return lnPIrwsat_wrapper_(this);
  }

  struct lnPIrwnmx_wrapper_ {
    explicit lnPIrwnmx_wrapper_(CriteriaWLTMMC* this_)
    : this_(this_) {}
    double operator ( )(const double & value) {
      return this_->lnPIrwnmx_(value);
    }
    CriteriaWLTMMC* this_;
  };

  lnPIrwnmx_wrapper_ lnPIrwnmx_wrap_(void) {
    return lnPIrwnmx_wrapper_(this);
  }

  // clone design pattern
  virtual shared_ptr<Criteria> cloneImpl_() const;
};

/// Factory method.
shared_ptr<CriteriaWLTMMC> makeCriteriaWLTMMC(const double beta,
  const double activ, const char* mType,
  const double mMin, const double mMax, const int nBin);

/// Factory method.
shared_ptr<CriteriaWLTMMC> makeCriteriaWLTMMC(const double beta,
  const double activ, const char* mType,
  const int nMin, const int nMax);

/// Factory method.
shared_ptr<CriteriaWLTMMC> makeCriteriaWLTMMC(const argtype &args);

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // CRITERIA_WLTMMC_H_
