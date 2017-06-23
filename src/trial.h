/**
 * \file
 *
 * \brief trial moves for Monte Carlo
 *
 */

#ifndef TRIAL_H_
#define TRIAL_H_

#include "./base_random.h"
#include "./space.h"
#include "./pair.h"
#include "./criteria.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

class Trial : public BaseRandom {
 public:
  Trial();
  Trial(Space *space, Pair* pair, Criteria* criteria);
  Trial(Space *space, Pair* pair, Criteria* criteria, const char* fileName);
  virtual ~Trial() {}
  virtual Trial* clone(Space* space, Pair* pair, Criteria* criteria) const = 0;
  shared_ptr<Trial> cloneShrPtr(Space* space, Pair* pair, Criteria* criteria) {
    return cloneImpl(space, pair, criteria); }

  // default constructor
  void defaultConstruction();

  // reset object pointers
  void reconstruct(Space* space, Pair *pair, Criteria *criteria);

  /// write restart file
  virtual void writeRestart(const char* fileName)
    { writeRestartBase(fileName); }
  void writeRestartBase(const char* fileName);

  // replace or restore criteria pointer
  void replaceCriteria(Criteria *criteria);
  void restoreCriteria();

  /// attempt trial
  virtual void attempt1() = 0;
  void attempt();

  /// initialize trial move counters and statistics
  void zeroStat();

  void trialAccept();     //!< call when accepting a trial
  void trialReject();     //!< call when rejecting a trial

  /// set the number of first bead attempts
  void numFirstBeads(const int nf);

  /// compute rosenbluth weights of multiple first bead attempts
  //    and select one of the trials
  // flag==0, old configuration (no moves, ins or dels)
  // flag==1, new configuration, same number of particles
  // flag==2, old configuration, preparing to delete (same as 0)
  // flag==3, just inserted particle
  double multiFirstBead(const int flag);

  /// initialize aggregation volume bias insertions or deletions
  // Aggregation-volume-bias Monte Carlo simulations of vapor-liquid
  // nucleation barriers for Lennard-Jonesium
  // J. Chem. Phys., Vol. 115, No. 23, 15 December 2001
  void initAVB(const double rAbove, const double rBelow);

  // trial move routines
  /// call when recording old configuration
  void trialMoveRecord();

  /// call when recording old configuration before collective move
  void trialMoveRecordAll(const int flag);

  /// call to decide whether to accept or reject trial move
  void trialMoveDecide(const double def, const double preFac);

  // tune parameters (e.g., based on acceptance)
  virtual void tuneParameters() {}
  void updateMaxMoveParam(const double percent, const double upperLimit,
    const double lowerLimit, const double targAcceptPer);
  double maxMoveParam;  //!< maximum parameter for move
  int maxMoveFlag;      //!< maxMoveParam is changing if == 1 (print to log)

  /// factory method
  shared_ptr<Trial> makeTrial(Space* space, Pair* pair, Criteria* criteria,
                              const char* fileName);

  /// return status of trial
  virtual string printStat(const bool header = false);

  // functions for read-only access of private data-members
  long long attempted() const { return attempted_; }   //!< number of attempts
  long long accepted() const { return accepted_; }
  double acceptPer() const;   //!< return trial acceptance percentage
  double de() const { return de_; }    //!< change in energy from last trial
  double deTot() const { return deTot_; }
  int nf() const { return nf_; }
  bool avbOn() const { return avbOn_; }
  Criteria* criteria() const { return criteria_; }
  double rAbove() const { return rAbove_; }
  double rBelow() const { return rBelow_; }

 protected:
  Space *space_;
  Pair *pair_;
  Criteria *criteria_;        //!< acceptance criteria
  string trialType_;          //!< type of trial
  Criteria *criteriaOld_;     //!< old acceptance criteria
  double de_;                 //!< change in energy due to trial move
  double peOld_;
  double preFac_;             //!< prefactor for acceptance criteria
  int reject_;                //!< outright reject move
  long double lnpMet_;        //!< metropolis probability
  double deTot_;              //!< change in energy due to multiple trial moves
  long long accepted_;    //!< counts number of accepted trials
  long long attempted_;   //!< counts number of attempted trials
  vector<int> mpart_;    //!< particles which are changed in the trial

  // multiple first bead variables
  double def_;                //!< change in energy due to first bead
  int nf_;                    //!< number of first bead trials
  vector<double> en_;    //!< energy of first bead trials
  vector<double> w_;     //!< Rosenbluth factor of of first bead trials
  vector<double> cpdf_;  //!< cumulative probability distribution function
                         //!< to select from Rosenbluth factors

  // avb variables
  double rAbove_;        //!< upper limit of bond
  double rBelow_;        //!< lower limit of bond
  double vIn_;           //!< aggregation volume
  bool avbOn_;           //!< aggregation volume bias flag
  std::string region_;   //!< selected avb region to move as bonded or nonbonded
  vector<int> tmpart_;   //!< target for avb move

  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl(Space* space, Pair *pair,
                                      Criteria *criteria) const  = 0;
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // TRIAL_H_

