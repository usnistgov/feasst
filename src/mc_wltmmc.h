/**
 * \file
 *
 * \brief attempts monte carlo trials and computes appropriate quantities
 *
 */

#ifndef WLTMMC_H_
#define WLTMMC_H_

#include "./mc.h"

namespace feasst {

class WLTMMC : public MC {
 public:
  WLTMMC(Space* space, Pair* pair, CriteriaWLTMMC* criteria);
  explicit WLTMMC(const char* fileName);
  ~WLTMMC();
  virtual WLTMMC* clone() const;
  shared_ptr<WLTMMC> cloneShrPtr() const
    { return std::static_pointer_cast<WLTMMC, MC>(cloneImpl()); }
  virtual void reconstruct();
  void defaultConstruction();

  /// write restart file
  void writeRestart(const char* fileName);

  /// run for a number of sweeps
  void runNumSweeps(const int nSweeps, const long long nprMax);
  void runNumSweepsExec(const int t, const int nSweeps,
                        vector<shared_ptr<WLTMMC> > &clones);
  void runNumSweepsRestart(const int nSweeps, const char* fileName);

  /// initialize collection matrix file name
  void initColMat(const char* fileName, const int nfreq)
    { colMatFileName_.assign(fileName); nFreqColMat_ = nfreq; };

  /// initialize OMP windows
  void initWindows(int flag) { if (flag == 1) { window_ = true; nExp_ = 1.5;
                  nOverlap_ = 2; nWindow_ = 1;} else { window_ = false; } }
  void initWindows(const double nExp, const int nOverlap)
    { window_ = true; nExp_ = nExp; nOverlap_ = nOverlap; nWindow_ = 1; }
  void initWindowsParaTemp(const double betaInc, const double lnzInc)
    { window_ = true; nWindow_ = 1; betaInc_ = betaInc; lnzInc_ = lnzInc; }

  /// this function is called after every trial attempt
  void afterAttempt();

  /// determine maximum number of particles for a given temperature
  //   and large activity
  int nMolMax(const long long npr, const double activ, const int nMolExtra);
  int nMolMax(const long long npr, const double activ)
    { return nMolMax(npr, activ, 0); }

  /// resize the nmol window in WLTMMC criteria based on lnPI or lnPIrw
  void nMolResizeWindow(const double liquidDrop, const int round);

  /// seek particle number which is in the range of WLTMMC
  void nMolSeekInRange(const int nMin, const int nMax);
  void nMolSeekInRange() { nMolSeekInRange(-1, -1); }

  /// append to all fileNames
  void appendFileNames(const char* chars) { colMatFileName_.append(chars);
    GRFileName_.append(chars); MC::appendFileNames(chars); }

  /// check that criteria of all trials match
  int checkTrialCriteria();

  /// print saturation summary
  vector<double> printSat();

  /// initialize overlapping processors
  void initOverlaps(const int t, vector<shared_ptr<WLTMMC> > &clones);

  /// initialize density threshold for Configurational Bias,
  //   only performed when nMolMax/V > thres
  void initConfigBiasDensThres(const double thres)
    { densThresConfigBias_ = thres; }
  void initNMolSeekTarget(const int target) { nMolSeekTarget_ = target; }

  /// add configuration swap trial
  #if defined (MPI_H_) || (_OPENMP)
    void confSwapTrial() { MC::confSwapTrial();
      trialConfSwapVec_.back()->initMType(c_->mType().c_str()); }
  #endif  // MPI_H_ || _OPENMP

  /// initialize GR file name
  void initGR(const char* fileName, const int nfreq, const double dr,
              const int iType, const int jType)
    { GRFileName_.assign(fileName); nFreqGR_ = nfreq;
      grt_.push_back(make_shared<Histogram>(dr, iType, jType));
      gr_.resize(grt_.size()); }
  void initGR(const char* fileName, const int nfreq, const double dr)
    { initGR(fileName, nfreq, dr, 0, 0); }

  /// initialize production run
  void wlFlatProduction(const int wlFlatProd) { wlFlatProd_ = wlFlatProd; }
  void wlFlatTerminate(const int wlFlat) { wlFlatTerm_ = wlFlat; }

  /// exposed pointers for easy access with clones
  CriteriaWLTMMC* c() { return c_; }
  int nWindows() const { return nWindow_; }

 protected:
  CriteriaWLTMMC* c_;           //!< wltmmc acceptance criteria
  std::string colMatFileName_;  //!< collection matrix file name
  int nFreqColMat_;  //!< frequency to update and print collection matrix

  /// production begins when wlFlat reaches this many (never if == -1)
  int wlFlatProd_;

  /// production ends when wlFlat reaches this many (never if == -1)
  int wlFlatTerm_;

  // OMP parallel windowing variables
  bool window_;   //!< windowing on or off
  int nWindow_;   //!< number of windows (default 1, or OMP_NUM_THREADS)
  double nExp_;      //!< window exponential scaling
  int nOverlap_;  //!< window overlap
  double betaInc_;  //!< parallel tempering if != 0
  double lnzInc_;  //!< parallel tempering if != 0

  // configurational bias flags
  double densThresConfigBias_;  //!< set configurational bias density threshold
  int nMolSeekTarget_;          //!< set number of molecules

  // radial distribution function (GR)
  string GRFileName_;         //!< GR file name
  int nFreqGR_;               //!< frequency to print GR

  /// canonical ensemble radial distribution functions for
  //   various pairs of particle types (e.g. gr_[type][nBins])
  vector<vector<shared_ptr<Histogram> > > gr_;

  /// gr templated used in initialization
  vector<shared_ptr<Histogram> > grt_;

  // clone design pattern
  virtual shared_ptr<MC> cloneImpl() const;
};

}  // namespace feasst

#endif  // WLTMMC_H_

