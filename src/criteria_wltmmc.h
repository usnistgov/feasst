/**
 * \file
 *
 * \brief acceptance criteria for Wang-Landau (WL), Transition Matrix (TM)
 * monte carlo (MC)
 *
 */

#ifndef CRITERIAWLTMMC_H_
#define CRITERIAWLTMMC_H_

#include "./criteria.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

class CriteriaWLTMMC : public Criteria {
 public:
  CriteriaWLTMMC(const double beta, const double activ, const char* mType,
                 const double mMin, const double mMax, const int nBin);
  explicit CriteriaWLTMMC(const char* fileName);
  ~CriteriaWLTMMC() {}
  CriteriaWLTMMC* clone() const;
  shared_ptr<CriteriaWLTMMC> cloneShrPtr() const;

  /// defaults in constructor
  void defaultConstruction();

  /// write restart file
  void writeRestart(const char* fileName);

  /// initialize the macrostate bins
  void initBins(const double mMax, const double mMin, const int nBin);

  /// acceptance criteria for trial moves
  int accept(const double pMet, const double peNew, const char* moveType,
             const int reject);

  /// store macrostate variables of old configuration
  void store(const Space* space, Pair* pair);

  /// given macrostate m, return bin number, or vice versa
  double bin2m(const int bin) const { return mMin_ + (bin + 0.5)*mBin_; }
  int bin(const double m) const { return feasstRound((m - bin2m(0))/mBin_); }

  /// set the modification factor
  void setg(const double g) { g_ = g; }

  /// if flatness criteria is met, reset histogram and reduce update factor
  int flatCheck();

  /// decide whether to use WL acceptance, TMMC acceptence, and whether to
  //  update collection matrix
  void collectInit() {collect_ = true; }
  void collectInit(const double lnfCollect) {lnfCollect_ = lnfCollect; }
  void tmmcInit(const double lnfTMMC) {lnfTMMC_ = lnfTMMC; }
  void tmmcInit();

  /// update histograms and collection matrix for visited macrostate
  void mUpdate(const int mOldBin, const int mNewBin, const double pmet,
               const int acceptFlag);

  /// update macrostate probability distribution with collection matrix
  void lnPIupdate();    // single simulation, local to this class

  /// sum multiple collection matrices
  void lnPIupdate(const vector<std::shared_ptr<CriteriaWLTMMC> > &c);

  /// convert collection matrix to lnPI
  void c2lnPI(const vector<vector<long double> > &col,
              vector<long double> *lnpiPtr);

  /// area under lnPI
  template<class T>
  double lnPIarea(const vector<T> &lnPI) const {
    double sum = 0;
    for (unsigned int i = 0; i < lnPI.size(); ++i) sum += exp(lnPI[i]);
    return sum;
  }
  double lnPIarea() const { return lnPIarea(lnPI_); }
  double lnPIrwarea() const { return lnPIarea(lnPIrw_); }

  /// normalize lnPI such that Sum(Pi) = 1
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
  void lnPInorm() { lnPInorm(&lnPI_); }

  /// average macrostate from lnPI
  template<class T>
  double lnPIaverage(const vector<T> &lnPI) {
    double av = 0;
    for (unsigned int m = 0; m < lnPI.size(); ++m) {
      av += bin2m(m) * exp(lnPI[m]);
    }
    return av/lnPIarea(lnPI);
  }
  double lnPIaverage() { return lnPIaverage(lnPI_); }
  double lnPIrwaverage() { return lnPIaverage(lnPIrw_); }

  /// set phase boundary
  void setPhaseBoundary(const int index) { phaseBoundary_ = index; }

  /// return phase boundaries, e.g. lnPI minima that are not first and last
  //  points
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
  vector<int> lnPIphaseBoundary() { return lnPIphaseBoundary(lnPI_); }
  vector<int> lnPIrwphaseBoundary() { return lnPIphaseBoundary(lnPIrw_); }

  /// number of phases
  template<class T>
  int lnPInumPhases(const vector<T> &lnPI) {
    vector<int> min = lnPIphaseBoundary(lnPI);
    return 1 + static_cast<int>(min.size());
  }
  int lnPInumPhases() { return lnPInumPhases(lnPI_); }
  int lnPIrwnumPhases() { return lnPInumPhases(lnPIrw_); }

  /// split phases among vector of criteria
  vector<CriteriaWLTMMC> phaseSplit(const vector<long double> &lnPI);

  /// split lnPI among phases
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


  /// print or read collection matrix and lnPI
  void printCollectMat(const char* fileName);

  /// paste together windows of collection matrices to print final version
  void printCollectMat(const char* fileName, const vector<CriteriaWLTMMC*> c);
  void printCollectMat(const char* fileName,
    const vector<shared_ptr<CriteriaWLTMMC> > c) {
    return printCollectMat(fileName, shrPtr2Raw(c));
  }

  /// read collection matrix
  void readCollectMat(const char* fileName);
  void readlnPIEner(const char* fileName);
  void readlnPIEnerCol(const char* fileName);

  /// zero all statistics and accumulators
  void zeroStat();

  /// reweight lnPI to different value of activity
  void lnPIrw(const double activrw);

  /// find saturation by reweighting and minimizing squared difference of peak
  //  heights
  double lnPIrwsat(const double activrw);

  /// find peak location by reweighting and minimizing squared different in
  //  peak location
  double lnPIrwnmx(const double lnactivrw);

  /// reweight to saturation conditions my mimiizing lnPIrw2sat
  void findSat();
  void printRWinit() { printRW_ = true; }

  /// reweight to obtain peak maxima at given nmol
  void findPeak(const double nMolPeak, const double lnactivGuess);

  /// reweight to find difference in peaks
  void peakDiff(const int iPeak, const int jPeak);

  /// return new nMolMax to resize the nmol window in WLTMMC criteria based on
  //  lnPI or lnPIrw
  int nMolResizeWindow(const double liquidDrop, const int round);

  /// pressure from lnPI
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

  /// obtain pressure isotherm from lnPI
  void lnPIpressureIso(const double vol);

  /// obtain grand canonical ensemble average from canonical ensemble average
  vector<double> lnPIgc2can(vector<double> data);
  void lnPIgc2can(const char* fileNameIn, const char* fileNameOut);

  /// data from lnPI
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

  /// obtain pressure isotherm from lnPI
  void lnPIenergyIso();

  /// replace lnPI
  void lnPIreplace(const vector<long double> &lnPI);

  /// splice windows together
  void spliceWindows(const vector<CriteriaWLTMMC*> c);
  void spliceWindows(const vector<shared_ptr<CriteriaWLTMMC> > c)
    { return spliceWindows(shrPtr2Raw(c)); }

  /// minimum number of sweeps in all windows
  int minNSweep(const vector<CriteriaWLTMMC*> c);
  int minNSweep(const vector<shared_ptr<CriteriaWLTMMC> > c)
    { return minNSweep(shrPtr2Raw(c)); }

  /// minimum number of wlFLat in all windows
  int minNwlFlat(const vector<CriteriaWLTMMC*> c);
  int minNwlFlat(const vector<shared_ptr<CriteriaWLTMMC> > c)
    { return minNwlFlat(shrPtr2Raw(c)); }

  /// prefill collection matrix with a constant
  void prefilColMat(const long double constant);

  /// return heat capacity for each bin
  vector<double> heatCapacity();
  vector<double> fluct();

  /// functor wrappers to pass to numerical recipe minimization algorithms
  struct lnPIrwsatwrapper {
    explicit lnPIrwsatwrapper(CriteriaWLTMMC* this_)
    : this_(this_) {}
    double operator ( )(const double & value) {
      return this_->lnPIrwsat(value);
    }
    CriteriaWLTMMC* this_;
  };

  lnPIrwsatwrapper lnPIrwsatwrap() {
    return lnPIrwsatwrapper(this);
  }

  struct lnPIrwnmxwrapper {
    explicit lnPIrwnmxwrapper(CriteriaWLTMMC* this_)
    : this_(this_) {}
    double operator ( )(const double & value) {
      return this_->lnPIrwnmx(value);
    }
    CriteriaWLTMMC* this_;
  };

  lnPIrwnmxwrapper lnPIrwnmxwrap(void) {
    return lnPIrwnmxwrapper(this);
  }

  /// Return current macrostate updated every call to accept() or every trial.
  int iMacro() const { return bin(mNew_); }

  /// read-only functions for protected variables
  string mType() const { return mType_; }
  double lnf() const { return lnf_; }
  vector<long double> lnPI() { return lnPI_; }
  vector<long double> lnPIrw() const { return lnPIrw_; }
  vector<long double> lnPIaggre() const { return lnPIaggre_; }
  vector<vector<long double> > C() const { return C_; }
  double lastbin2m() const { return bin2m(nBin_ - 1); }
  double mBin() const { return mBin_; }
  double mMax() const { return mMax_; }
  double mMin() const { return mMin_; }
  double mOld() const { return mOld_; }
  double mNew() const { return mNew_; }
  int nBin() const { return nBin_; }
  double nTunnels() const { return nTunnels_; }
  double nSweep() const { return nSweep_; }
  double activrw() const { return activrw_; }
  bool printRW() const { return printRW_; }
  bool collect() const { return collect_; }
  bool tmmc() const { return tmmc_; }
  vector<double> lnpi2pressure() const { return lnpi2pressure_; }
  vector<Accumulator> pe() const { return pe_; }
  int wlFlat() const { return wlFlat_; }
  double lnfCollect() const { return lnfCollect_; }
  vector<double> peMUVT() { lnPIenergyIso(); return peMUVT_; }

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

  // clone design pattern
  virtual shared_ptr<Criteria> cloneImpl_() const;
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // CRITERIAWLTMMC_H_

