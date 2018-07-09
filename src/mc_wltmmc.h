/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef WLTMMC_H_
#define WLTMMC_H_

#include <memory>
#include <string>
#include <vector>
#include "./mc.h"

namespace feasst {

/**
 * Attempts Monte Carlo trials with flat-histogram methods.
 */
class WLTMMC : public MC {
 public:
  /// Constructor
  WLTMMC(Space* space, Pair* pair, CriteriaWLTMMC* criteria);

  /// Constructor
  WLTMMC(shared_ptr<Space> space, shared_ptr<Pair> pair,
    shared_ptr<CriteriaWLTMMC> criteria)
    : WLTMMC(space.get(), pair.get(), criteria.get()) {}

  /// Constructor
  WLTMMC(shared_ptr<Pair> pair, shared_ptr<CriteriaWLTMMC> criteria)
    : WLTMMC(pair->space(), pair.get(), criteria.get()) {}

  /// Initialize collection matrix file name.
  void initColMat(const char* fileName, const int nfreq)
    { colMatFileName_.assign(fileName); nFreqColMat_ = nfreq; };

  /// Run for a number of sweeps.
  void runNumSweeps(const int nSweeps, const long long nprMax = -1);

  /// Initialize production run.
  void wlFlatProduction(const int wlFlatProd) { wlFlatProd_ = wlFlatProd; }

  /// Set this to 1 if you want to use runNumSweeps(wlFlat) where you input
  /// wlFlats instead of sweeps.
  void wlFlatTerminate(const int wlFlat) { wlFlatTerm_ = wlFlat; }

  /** For renaming files for individual processors, set the appended name.
   *  For example, the default "p" would rename files "file" -> "filep0" */
  void setProcessorFileDescription(const char* append = "_core") {
    procFileAppend_.assign(append);
  }

  /// Run number of sweeps from restart file which is the clone for a spawn
  /// of parallel runs.
  void runNumSweepsRestart(const int nSweeps, const char* fileName);

  /// Initialize OMP windows.
  /// Note that the number of OMP processors is set by the bash environment.
  /// For example, "export OMP_NUM_THREADS=4".
  void initWindows(int flag   //!< use OMP if flag == 1
    ) { if (flag == 1) { window_ = true; nExp_ = 1.5;
    nOverlap_ = 0; nWindow_ = 1;} else { window_ = false; } }

  /// Initialize OMP windows with nExp determining spacing and nOverlap extra
  /// macrostates with processor overlaps.
  void initWindows(const double nExp, const int nOverlap)
    { window_ = true; nExp_ = nExp; nOverlap_ = nOverlap; nWindow_ = 1; }

  /// Initialize OMP windows for parallel tempering.
  void initWindowsParaTemp(const double betaInc, const double lnzInc)
    { window_ = true; nWindow_ = 1; betaInc_ = betaInc; lnzInc_ = lnzInc; }

  // determine maximum number of particles for a given temperature
  //   and large activity
  int nMolMax(const long long npr, const double activ, const int nMolExtra);
  int nMolMax(const long long npr, const double activ)
    { return nMolMax(npr, activ, 0); }

  /// Seek particle number which is in the range of WLTMMC.
  void nMolSeekInRange(const int nMin, const int nMax);
  void nMolSeekInRange() { nMolSeekInRange(-1, -1); }

  /// Append chars to all file names.
  void appendFileNames(const char* chars) { colMatFileName_.append(chars);
    MC::appendFileNames(chars); }

  /// Check that criteria of all trials match.
  int checkTrialCriteria();

// HWH mins
//  /// Print saturation summary.
//  vector<double> printSat();

  /// Initialize overlapping processors for OMP configuration swap trials.
  void initOverlaps(const int t, vector<shared_ptr<WLTMMC> > *clones);

  /// Initialize density threshold such that Configurational Bias is
  /// only performed when nMolMax/V > thres.
  void initConfigBiasDensThres(const double thres)
    { densThresConfigBias_ = thres; }
  void initNMolSeekTarget(const int target) { nMolSeekTarget_ = target; }

  #if defined (MPI_H_) || (_OPENMP)
    /// Add configuration swap trial.
    void confSwapTrial() { MC::confSwapTrial();
      trialConfSwapVec_.back()->initMType(c_->mType().c_str()); }
  #endif  // MPI_H_ || _OPENMP

  /// Return pointer to WLTMMC acceptance criteria class
  CriteriaWLTMMC* c() { return c_; }

  /// Return the number of windows for OMP parallelization.
  int nWindows() const { return nWindow_; }

  /// Write restart file.
  void writeRestart(const char* fileName);

  /// Construct from restart file.
  explicit WLTMMC(const char* fileName);

  ~WLTMMC();
  virtual WLTMMC* clone() const;
  shared_ptr<WLTMMC> cloneShrPtr() const
    { return std::static_pointer_cast<WLTMMC, MC>(cloneImpl()); }
  virtual void reconstruct();

 protected:
  CriteriaWLTMMC* c_;           //!< wltmmc acceptance criteria
  std::string colMatFileName_;  //!< collection matrix file name
  std::string procFileAppend_;
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

  /// Execute running number of sweeps.
  void runNumSweepsExec_(const int t, const int nSweeps,
                        vector<shared_ptr<WLTMMC> > *clones);

  /// this function is called after every trial attempt
  void afterAttempt_();

  void defaultConstruction_();

  // clone design pattern
  virtual shared_ptr<MC> cloneImpl() const;
};

}  // namespace feasst

#endif  // WLTMMC_H_
