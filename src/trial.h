#ifndef TRIAL_H_
#define TRIAL_H_

#include <memory>
#include <string>
#include <vector>
#include "./base_random.h"
#include "./space.h"
#include "./pair.h"
#include "./criteria.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Trial move for Monte Carlo simulations.
 */
class Trial : public BaseRandom {
 public:
  /// Constructor.
  Trial(Space *space, Pair* pair, Criteria* criteria);

  // Reset object pointers.
  void reconstruct(Space* space, Pair *pair, Criteria *criteria);

  /// This constructor is not often used, but its purpose is to initialize trial
  /// for interface before using reconstruct to set object pointers.
  Trial();

  /// Replace or restore criteria pointer.
  void replaceCriteria(Criteria *criteria);
  void restoreCriteria();

  /// Attempt trial.
  void attempt();

  /// Initialize trial move counters and statistics.
  void zeroStat();

  /// Set the number of first bead attempts, nf.
  void numFirstBeads(const int nf);

  /**
   * Initialize aggregation volume bias region.
   * Aggregation-volume-bias Monte Carlo simulations of vapor-liquid
   * nucleation barriers for Lennard-Jonesium
   * J. Chem. Phys., Vol. 115, No. 23, 15 December 2001
   * @param rAbove upper limit of bond
   * @param rBelow lower limit of bond
   */
  void initAVB(const double rAbove, const double rBelow);

  double maxMoveParam;  //!< maximum parameter for move
  int maxMoveFlag;      //!< maxMoveParam is changing if == 1 (print to log)

  // Tune parameters based on percentage acceptance.
  virtual void tuneParameters() {}

  // Note: this can be moved to a class only involves TrialAdd, TrialDelete
  /**
   * Confine trials to particular values of a given dimension.
   * Note that these confines apply to the center of mass.
   * @param upper upper limit of insertions
   * @param lower lower limit of insertions
   * @param dimension dimension for bounds above (e.g., z is "2")
   */
  void confine(const double upper, const double lower, const int dimension);

  /// Indicate that the trial is being performed to reach a certain number of
  /// particles, nPartTarget. Therefore, detailed balance may not need to be
  /// obeyed in order to speed up nMolSeek.
  /// If nPartTarget == -1, then there is no seeking.
  void initializeNMolSeek(const int nPartTarget = -1) {
    nPartTarget_ = nPartTarget; }

  /// Factory method.
  shared_ptr<Trial> makeTrial(Space* space, Pair* pair, Criteria* criteria,
                              const char* fileName);

  /// Return status of trial
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

  /// Write restart file.
  virtual void writeRestart(const char* fileName)
    { writeRestartBase(fileName); }
  void writeRestartBase(const char* fileName);

  /// Construct from restart file.
  Trial(Space *space, Pair* pair, Criteria* criteria, const char* fileName);

  virtual ~Trial() {}
  virtual Trial* clone(Space* space, Pair* pair, Criteria* criteria) const = 0;
  shared_ptr<Trial> cloneShrPtr(Space* space, Pair* pair, Criteria* criteria) {
    return cloneImpl(space, pair, criteria); }

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

  virtual void attempt1_() = 0;  // this is the implementation of a trial.

  void trialAccept_();     //!< call when accepting a trial
  void trialReject_();     //!< call when rejecting a trial

  /// call when recording old configuration
  void trialMoveRecord_();

  /// call when recording old configuration before collective move
  void trialMoveRecordAll_(const int flag);

  /// call to decide whether to accept or reject trial move
  /// @param def energy of first bead
  /// @param preFac acceptance criteria prefactor
  void trialMoveDecide_(const double def, const double preFac);

  /// Update max move parameters.
  void updateMaxMoveParam_(const double percent, const double upperLimit,
    const double lowerLimit, const double targAcceptPer);

  // multiple first bead variables
  double def_;                //!< change in energy due to first bead
  int nf_;                    //!< number of first bead trials
  vector<double> en_;    //!< energy of first bead trials
  vector<double> w_;     //!< Rosenbluth factor of of first bead trials
  vector<double> cpdf_;  //!< cumulative probability distribution function
                         //!< to select from Rosenbluth factors

  /// compute rosenbluth weights of multiple first bead attempts
  ///    and select one of the trials.
  /// flag==0, old configuration (no moves, ins or dels)
  /// flag==1, new configuration, same number of particles
  /// flag==2, old configuration, preparing to delete (same as 0)
  /// flag==3, just inserted particle
  double multiFirstBead_(const int flag);

  // avb variables
  double rAbove_;        //!< upper limit of bond
  double rBelow_;        //!< lower limit of bond
  double vIn_;           //!< aggregation volume
  bool avbOn_;           //!< aggregation volume bias flag
  std::string region_;   //!< selected avb region to move as bonded or nonbonded
  vector<int> tmpart_;   //!< target for avb move

  // confinement variables
  int confineFlag_ = 0;   //!< set to "1" if confinement is on
  double confineUpper_;
  double confineLower_;
  int confineDim_;

  int nPartTarget_;   //!< target number of particles

  // default constructor
  void defaultConstruction_();

  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl(Space* space, Pair *pair,
                                      Criteria *criteria) const  = 0;
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // TRIAL_H_

