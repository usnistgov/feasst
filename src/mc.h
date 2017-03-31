/**
 * \file
 *
 * \brief attempts monte carlo trials and analyzes quantities
 *
 */

#ifndef MC_H_
#define MC_H_

#include "./criteria_metropolis.h"
#include "./criteria_wltmmc.h"
#include "./trial.h"
#ifdef MPI_H_
  #include "./trial_confswap_txt.h"
#endif  // MPI_H_
#ifdef OMP_H_
  #include "./trial_confswap_omp.h"
#endif  // OMP_H_
#include "./accumulator.h"
#include "./analyze.h"

class MC : public BaseAll {
 public:
  MC(Space *space, Pair *pair, Criteria *criteria);
  explicit MC(const char* fileName);
  virtual ~MC();
  virtual MC* clone() const;
  shared_ptr<MC> cloneShrPtr() const { return cloneImpl(); }
  shared_ptr<MC> cloneShallowShrPtr() const { return cloneShallowImpl(); }
  virtual void reconstruct();
  void defaultConstruction();

  /// write restart file
  virtual void writeRestart(const char* fileName);

  /// attempt trial move according to weights
  void attemptTrial();

  /// add configuration swap trial
  virtual void confSwapTrial();

  /// initialize trial
  void initTrial(Trial* trial);
  void initTrial(shared_ptr<Trial> trial);

  /// remove trial
  void removeTrial(const int iTrial);

  /// zero statistics of all trials
  void zeroStat();

  double weight;                  //!< weight of trials, default 1

  /// function to quickly seek nMol particles as initial configuration
  void nMolSeek(const int n, const char* molType, long long maxAttempts);
  void nMolSeek(const int n, long long maxAttempts)
    { nMolSeek(n, "", maxAttempts); }
  void nMolSeek(const int n, const char* molType)
    { nMolSeek(n, molType, 1e12); }

  /// determine maximum number of particles for a given temperature
  //   and large activity
  virtual int nMolMax(const long long npr, const double activ,
    const int nMolExtra);
  virtual int nMolMax(const long long npr, const double activ)
    { return nMolMax(npr, activ, 0); }

  /// set the number of trials in a production simulation
  void setNumTrials(const long long npr) { npr_ = npr; }
  void runNumTrials(const long long npr) { setNumTrials(npr); run(); }

  /// run a production level simulation
  virtual void run();

  /// print functions
  void printStat();     //!< print status of all trials to log
  double pePerMol();     //!< print potential energy per molecule

  /// this function is called after every trial attempt
  //   derived classes must call afterAttemptBase()
  virtual void afterAttempt() { afterAttemptBase(); }
  void afterAttemptBase();

  /// turn on neigh list for avb trials,
  //   or check that it is on with matching region
  void neighAVBInit(const double rAbove, const double rBelow);

  /// initialize log file name
  void initLog(const char* fileName, const long long nfreq)
    { logFileName_.assign(fileName); nFreqLog_ = nfreq; }

  /// initialize movie file name
  void initMovie(const char* fileName, const int nfreq)
    { movieFileName_.assign(fileName); nFreqMovie_ = nfreq; }

  /// initialize XTC file name
  void initXTC(const char* fileName, const int nfreq)
    { XTCFileName_.assign(fileName); nFreqXTC_ = nfreq; }

  /// initialize Analyzer
  void initAnalyze(Analyze* analyze) { analyzeVec_.push_back(analyze); }

  /// append to all fileNames
  virtual void appendFileNames(const char* chars);
  virtual void appendProductionFileNames(const char* chars);

  /// initialize freqeuncy to check energy
  void setNFreqCheckE(const double nfreq, const double tol)
    { nFreqCheckE_ = nfreq; checkEtol_ = tol; }

  /// initialize frequency to tune trial transform parameters
  void setNFreqTune(const double nfreq) { nFreqTune_ = nfreq; }

  /// initialize restart file name
  void initRestart(const char* fileName, const int nfreq)
    { rstFileBaseName_.assign(fileName); rstFileName_.assign(fileName);
      nFreqRestart_ = nfreq; }

  /// check that criteria of all trials match
  virtual int checkTrialCriteria();

  /// compute second virial coefficient by Monte Carlo integration
  void b2(const double tol, double &b2v, double &b2er, double boxl);
  void b2(const double tol, double &b2v, double &b2er)
    { b2(tol, b2v, b2er, -1.); }
  vector<double> b2(const double tol)
    { vector<double> rtrn; double b2v, b2s; b2(tol, b2v, b2s);
      rtrn.push_back(b2v); rtrn.push_back(b2s); return rtrn; }

  /// compute the Boyle temperature, b2(T_Boyle)==0
  double boyle(const double tol);
  double boylemin(const double beta);

  /// functor wrappers to pass to numerical recipe minimization algorithms
  struct boyleminwrapper {
    explicit boyleminwrapper(MC* this_)
    : this_(this_) {}
    double operator ( )(const double & value) {
      return this_->boylemin(value);
    }
    MC* this_;
  };
  boyleminwrapper boyleminwrap(void) {
    return boyleminwrapper(this);
  }

  /// remove all configurational bias trials
  void removeConfigBias();

  /// update cumulative probability of trials
  void updateCumulativeProb();

  /// replace and restore pointer to criteria
  void replaceCriteria(Criteria *criteria);
  void restoreCriteria();

  /// initialize production run
  void initProduction();

  /// remove ownership of pointers
  void removeOwnership()
    { spaceOwned_ = false; pairOwned_ = false; criteriaOwned_ = false; }

  /// tune move parameters
  void tuneTrialParameters();

  /// read only access to protected variables
  int nTrials() const { return int(trialVec_.size()); }
  vector<shared_ptr<Trial> > trialVec() const { return trialVec_; }
  vector<double> trialWeight() const { return trialWeight_; }
  vector<double> trialCumulativeProb() const { return trialCumulativeProb_; }
  #ifdef MPI_H_
    TrialConfSwapMPI* trialConfSwap(const int i)
      { return trialConfSwapVec_[i].get(); }
  #endif  // MPI_H_
  #ifdef OMP_H_
    TrialConfSwapOMP* trialConfSwap(const int i)
      { return trialConfSwapVec_[i].get(); }
  #endif  // OMP_H_
  double enAv() const { return peStat_.average(); }
  double enStdev() const { return peStat_.stdev(); }
  double nMolAv() const { return nMolStat_.average(); }
  double nMolStdev() const { return nMolStat_.stdev(); }
  int id() const { return space_->id(); }
  long long nAttempts() const { return nAttempts_; }
  Criteria* criteria() const { return criteria_; }
  Space* space() const { return space_; }
  Pair* pair() const { return pair_; }
  string rstFileName() const { return rstFileName_; }
  Accumulator peStat() const { return peStat_; }
  long long nFreqLog() const { return nFreqLog_; }
  int nFreqMovie() const { return nFreqMovie_; }
  vector<Analyze*> analyzeVec() const { return analyzeVec_; }
  bool spaceOwned() const { return spaceOwned_; }

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
  #ifdef OMP_H_
    /// vector of ConfSwap trials
    vector<shared_ptr<TrialConfSwapOMP> > trialConfSwapVec_;
  #endif  // OMP_H_
  Accumulator peStat_;        //!< accumulate potential energy statistics
  Accumulator nMolStat_;      //!< accumulate number of molecule statistics
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
  long long npr_;         //!< number of trials in simulaiton
  double checkEtol_;          //!< tolerance for energy check
  bool printPressure_;        //!< flag to turn on printing of pressure
  int production_;           //!< flag for production simulation
  double boyletol_;           //!< tolerance for boyle temperature

  // analyzers
  vector<Analyze*> analyzeVec_;   //!< vector of pointers to analyzers

  // error messaging
  void mout_(const char* messageType, std::ostream& message)
    {myOut(messageType, message, className_, verbose_);}

  // clone design pattern
  virtual shared_ptr<MC> cloneImpl() const;
  virtual shared_ptr<MC> cloneShallowImpl() const;
};

#endif  // MC_H_

