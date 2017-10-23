/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#ifndef MC_H_
#define MC_H_

#include "./criteria_metropolis.h"
#include "./criteria_wltmmc.h"
#include "./trial.h"
#ifdef MPI_H_
  #include "./trial_confswap_txt.h"
#endif  // MPI_H_
#ifdef _OPENMP
  #include "./trial_confswap_omp.h"
#endif  // _OPENMP
#include "./accumulator.h"
#include "./analyze.h"
#include "./base_random.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Attempts Monte Carlo trials and analyzes quantities.
 */

class MC : public BaseRandom {
 public:
  /// Constructor
  MC(Space *space, Pair *pair, Criteria *criteria);

  /// Weight for probability of selection of Monte Carlo trials.
  double weight = 1.;

  /// Add trial to MC with trial weight set to current weight.
  void initTrial(shared_ptr<Trial> trial);

  /// Attempt a random trial according to weights.
  void attemptTrial();

  /// Add configuration swap trial.
  /// HWH: Depreciate, but requires communication between processors for OMP.
  virtual void confSwapTrial();

  /// Remove trial with index in order of initialization.
  /// If iTrial is not provided, remove the last trial that was added.
  void removeTrial(int iTrial = -1);

  // Set the number of trials in a simulation.
  void setNumTrials(const long long npr) { npr_ = npr; }

  // run a production level simulation
  virtual void run();

  /// Run a simulation with "npr" trials.
  void runNumTrials(const long long npr) { setNumTrials(npr); run(); }

  /// Initialize production run.
  void initProduction();

  /** For renaming files for production, set the appended name.
   *  For example, the default "pr" would rename files "file" -> "filepr" */
  void setProductionFileDescription(const char* append = "_prod") {
    prodFileAppend_.assign(append);
  }

  /// Zero statistics of all mc variables, criteria, and all trials.
  void zeroStat();

  /// Attempt to seek nMol particles of type molType by inserting/deleting
  /// or particles.
  /// Note that detailed balance is not obeyed.
  void nMolSeek(const int nMol, const char* molType,
    /// maximum number of trial attempts
    long long maxAttempts,
    /// If pressure is set in criteria and seeking more particles, then begin
    /// expanding the volume by a factor of volumeExpansion, insert/delete
    /// particles, then perform a volume trial move at high pressure until
    /// original box size is obtained.
    const double volumeExpansion = 1.2);
  void nMolSeek(const int nMol, long long maxAttempts = 1e12)
    { nMolSeek(nMol, "", maxAttempts); }
  void nMolSeek(const int nMol, const char* molType)
    { nMolSeek(nMol, molType, 1e12); }

  /// initialize log file name and number of trials per print.
  void initLog(const char* fileName, const long long nfreq)
    { logFileName_.assign(fileName); nFreqLog_ = nfreq; }

  /// initialize movie file name and number of trials per print.
  void initMovie(const char* fileName, const int nfreq)
    { movieFileName_.assign(fileName); nFreqMovie_ = nfreq; }

  /// Initialize XTC file name and number of trials per print.
  void initXTC(const char* fileName, const int nfreq)
    { XTCFileName_.assign(fileName); nFreqXTC_ = nfreq; }

  /// Initialize freqeuncy to check running energy against the total energy
  /// recalculated every nfreq trials. Error if not within tolerance.
  void setNFreqCheckE(const double nfreq, const double tolerance)
    { nFreqCheckE_ = nfreq; checkEtol_ = tolerance; }

  /// Initialize frequency to tune trial parameters.
  /// Note that tuning does not obey detailed balance.
  void setNFreqTune(const double nfreq) { nFreqTune_ = nfreq; }

  /// Initialize restart file name and print every nfreq trials.
  void initRestart(const char* fileName, const int nfreq)
    { rstFileBaseName_.assign(fileName); rstFileName_.assign(fileName);
      nFreqRestart_ = nfreq; }

  /// Initialize Analyzer.
	void initAnalyze(shared_ptr<Analyze> analyze) {
    analyze->reconstruct(pair_); analyzeVec_.push_back(analyze);
  }

  // determine maximum number of particles for a given temperature
  //   and large activity
  virtual int nMolMax(const long long npr, const double activ,
    const int nMolExtra);
  virtual int nMolMax(const long long npr, const double activ)
    { return nMolMax(npr, activ, 0); }

  /// print functions
  void printStat(const std::string hash="");     //!< print status of all trials to log
  double pePerMol();     //!< print potential energy per molecule

  /// turn on neigh list for avb trials,
  //   or check that it is on with matching region
  void neighAVBInit(const double rAbove,  //!< upper bound of spherical shell
                    const double rBelow   //!< lower bound of spherical shell
                   );

  /// Append chars to all fileNames.
  virtual void appendFileNames(const char* chars);

  /// Append chars to all file names associated with production.
  virtual void appendProductionFileNames(const char* chars);

  /// Check that criteria of all trials are the same.
  virtual int checkTrialCriteria();

  /**
   * compute second virial coefficient by Monte Carlo integration
   *   B2(T)=-0.5*int(dr*f(r)); f(r)=meyer fn
   *   B2(T)=-0.5*int(dr)*(1/npr)*sum[f(r)]
   *   B2(T)=-0.5*volume*(1/npr)*sum[f(r_i)]
   */
  void b2(const double tol, double &b2v, double &b2er, double boxl);
  void b2(const double tol, double &b2v, double &b2er)
    { b2(tol, b2v, b2er, -1.); }
  vector<double> b2(const double tol)
    { vector<double> rtrn; double b2v, b2s; b2(tol, b2v, b2s);
      rtrn.push_back(b2v); rtrn.push_back(b2s); return rtrn; }

  /**
   * Compute second virial coefficient by Mayer sampling Monte Carlo.
   * https://doi.org/10.1103/PhysRevLett.92.220601
   */
  void b2mayer(
    double *b2v,              //!< return value of the second virial coefficient
    double *b2er,             //!< standard deviation of the mean
    Pair *pairRef,            //!< reference potential
    const double tol = 1e-4,  //!< terminate trials when tolerance reached
    double boxl = -1          //!< box length containing all nonzero energy
  );

  /// compute the Boyle temperature, \f$ B_2(T_{Boyle})=0\f$.
  double boyle(const double tol);

  /// Remove all configurational bias trials.
  void removeConfigBias();

  /*
   * Replace pointer to criteria for MC and all trials.
   *
   * HWH WARNING: TEMPORARY/LIMITED USE CASE
   *
   * If one intends for this criteria to be present for the long term,
   * this criteria must be replaced, or ownership updated,
   * or else you'll get a memory leak when the destructor is called.
   */
  void replaceCriteria(Criteria *criteria);

  // Restore pointer to criteria.
  void restoreCriteria();

  /// remove ownership of pointers
  void removeOwnership()
    { spaceOwned_ = false; pairOwned_ = false; criteriaOwned_ = false; }

  /// read only access to protected variables
  int nTrials() const { return int(trialVec_.size()); }
  vector<shared_ptr<Trial> > trialVec() const { return trialVec_; }
  vector<double> trialWeight() const { return trialWeight_; }
  vector<double> trialCumulativeProb() const { return trialCumulativeProb_; }
  #ifdef MPI_H_
    TrialConfSwapMPI* trialConfSwap(const int i)
      { return trialConfSwapVec_[i].get(); }
  #endif  // MPI_H_
  #ifdef _OPENMP
    TrialConfSwapOMP* trialConfSwap(const int i)
      { return trialConfSwapVec_[i].get(); }
  #endif  // _OPENMP
  double enAv() const { return peAccumulator_.average(); }
  double enStdev() const { return peAccumulator_.stdev(); }
  double nMolAv() const { return nMolAccumulator_.average(); }
  double nMolStdev() const { return nMolAccumulator_.stdev(); }
  int id() const { return space_->id(); }
  long long nAttempts() const { return nAttempts_; }
  Criteria* criteria() const { return criteria_; }
  Space* space() const { return space_; }
  Pair* pair() const { return pair_; }
  string rstFileName() const { return rstFileName_; }
  Accumulator peAccumulator() const { return peAccumulator_; }
  long long nFreqLog() const { return nFreqLog_; }
  int nFreqMovie() const { return nFreqMovie_; }
  vector <shared_ptr <Analyze> > analyzeVec() const { return analyzeVec_; }
  bool spaceOwned() const { return spaceOwned_; }
  int production() const { return production_; }

  /// Write restart file.
  virtual void writeRestart(const char* fileName);

  /// Construction by restart file.
  explicit MC(const char* fileName);
  virtual ~MC();
  virtual MC* clone() const;
  shared_ptr<MC> cloneShrPtr() const { return cloneImpl(); }
  shared_ptr<MC> cloneShallowShrPtr() const { return cloneShallowImpl(); }
  virtual void reconstruct();

 protected:
  Space *space_;            //!< spatial information
  Pair *pair_;              //!< pairwise interactions
  Criteria *criteria_;      //!< acceptance criteria
  Criteria *criteriaOld_;   //!< old acceptance criteria
  bool spaceOwned_;         //!< flag if delete necessary on destruction
  bool pairOwned_;          //!< flag if delete necessary on destruction
  bool criteriaOwned_;      //!< flag if delete necessary on destruction
  vector<shared_ptr<Trial> > trialVec_;  //!< vector of trials
  vector<double> trialWeight_;           //!< vector of trial weights
  vector<double> trialCumulativeProb_;   //!< cumulative probability of trials
  #ifdef MPI_H_
    /// vector of ConfSwap trials
    vector<shared_ptr<TrialConfSwapTXT> > trialConfSwapVec_;
  #endif  // MPI_H_
  #ifdef _OPENMP
    /// vector of ConfSwap trials
    vector<shared_ptr<TrialConfSwapOMP> > trialConfSwapVec_;
  #endif  // _OPENMP
  Accumulator peAccumulator_;        //!< accumulate potential energy statistics
  Accumulator nMolAccumulator_;      //!< accumulate number of molecule statistics
  double prSum_;              //!< sum of pressure of each state
  double prSum2_;             //!< sum of square of pressure of each state
  long long nAttempts_;   //!< number of attempted trials

  /// flag to print header in log file
  //  if 0, no header
  //  if 1, print header
  //  if 2, print header with a comment "#" on first line (for restarts)
  //  if -1, print line with "#" but not header
  int printLogHeader_;

  string logFileName_;        //!< log file name
  long long nFreqLog_;    //!< frequency to print to log
  string movieFileName_;      //!< movie file name
  int nFreqMovie_;            //!< frequency to print movie
  string XTCFileName_;        //!< XTC file name
  int nFreqXTC_;              //!< frequency to print XTC
  int nFreqCheckE_;           //!< frequency to check energy
  int nFreqTune_;             //!< frequency to tune translation parameters
  int nFreqRestart_;          //!< frequency to write restart file
  string rstFileName_;        //!< restart file name
  string rstFileBaseName_;    //!< restart file base name
  std::string prodFileAppend_;
  long long npr_;         //!< number of trials in simulaiton
  double checkEtol_;          //!< tolerance for energy check
  bool printPressure_;        //!< flag to turn on printing of pressure
  int production_;           //!< flag for production simulation
  double boyletol_;           //!< tolerance for boyle temperature

  // analyzers
  vector<shared_ptr<Analyze> > analyzeVec_;

  // unique hash for configurations
  std::string hash_;

  // virial coefficient
  void b2init_();

  /// update cumulative probability of trials
  void updateCumulativeProb_();

  /// tune move parameters
  void tuneTrialParameters_();

  /// Return squared B2 to compute the boyle temperature.
  double boylemin_(const double beta);

  // functor wrappers to pass to numerical recipe minimization algorithms
  struct boyleminwrapper_ {
    explicit boyleminwrapper_(MC* this_)
    : this_(this_) {}
    double operator ( )(const double & value) {
      return this_->boylemin_(value);
    }
    MC* this_;
  };
  boyleminwrapper_ boyleminwrap_(void) {
    return boyleminwrapper_(this);
  }

  /// this function is called after every trial attempt
  //   derived classes must call afterAttemptBase_()
  virtual void afterAttempt_() { afterAttemptBase_(); }
  void afterAttemptBase_();

  void defaultConstruction_();

  // clone design pattern
  virtual shared_ptr<MC> cloneImpl() const;
  virtual shared_ptr<MC> cloneShallowImpl() const;
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // MC_H_
