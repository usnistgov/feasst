/**
 * \file
 *
 * \brief trial moves for Monte Carlo
 *
 */

#ifndef TRIAL_MD_H_
#define TRIAL_MD_H_

#include "./trial.h"

class TrialMD : public Trial {
 public:
  TrialMD();
  TrialMD(Space *space, Pair *pair, Criteria *criteria);
  TrialMD(const char* fileName, Space *space, Pair *pair,
             Criteria *criteria);
  ~TrialMD() {}
  TrialMD* clone(Space* space, Pair* pair, Criteria* criteria) const {
    TrialMD* t = new TrialMD(*this);
    t->reconstruct(space, pair, criteria); return t;
  }
  shared_ptr<TrialMD> cloneShrPtr(
    Space* space, Pair* pair, Criteria* criteria) const {
    return(std::static_pointer_cast<TrialMD, Trial>(
      cloneImpl(space, pair, criteria)));
  }
  void defaultConstruction();
  void writeRestart(const char* fileName);

  /// attempt random translation
  void attempt1();

  /// initialize momentum
  void initMomentum();

  /// zero total momentum of the system
  void zeroTotalMomentum();

  /// kinetic energy
  double kineticEnergy();

  /// instantaneous temperature from particle velocities
  double temperature();

  /// rescale velocity according to the temperature
  void rescaleVelocity();
  int rescaleTemp;     //!< rescale temperature (1==yes)

  double timestep;       //!< integration time step
  int nFreqZeroMomentum;   //!< zero momentum every this many steps

  /// return status of trial
  string printStat(const bool header = false);
 
  /// update center of mass forces
  void updateFCOM();

  /// update half-step velocities
  void updateVelocityHalfStep();
 
  /// velocity verlet integration
  void integrateVelocityVerlet();

  // read-only access to private
  double fCOM(const int iMol, const int dim) const { 
    return fCOM_[space_->dimen()*iMol + dim];
  }

protected:
  vector<double> vel_;    //!< particle velocities
  vector<double> mass_;   //!< particle masses (assumed unity)
  vector<double> fCOM_;   //!< center of mass forces
  string integrator_;     //!< type of integrator

  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl(
    Space* space, Pair *pair, Criteria *criteria) const {
    shared_ptr<TrialMD> t = make_shared<TrialMD>(*this);
    t->reconstruct(space, pair, criteria);
    return t;
  }
};

#endif  // TRIAL_MD_H_

