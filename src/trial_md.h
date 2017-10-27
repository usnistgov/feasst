/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef TRIAL_MD_H_
#define TRIAL_MD_H_

#include <memory>
#include <string>
#include <vector>
#include "./trial.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Perform a molecular dynamics (MD) integration step.
 * Note that only translational integration is performed.
 * So it is not correct for anisotropic or multi-site particles.
 */
class TrialMD : public Trial {
 public:
  /// Constructor
  TrialMD(Pair *pair, Criteria *criteria);

  /// This constructor is not often used, but its purpose is to initialize trial
  /// for interface before using reconstruct to set object pointers.
  TrialMD();

  double timestep;       //!< integration time step
  int nFreqZeroMomentum;   //!< zero momentum every this many steps
  int rescaleTemp;     //!< rescale temperature (1==yes)

  /**
   * Initialize momentum (the lazy way with plain rescaling).
   * Maxwell-Boltzmann is expected to rapidly estability after a few steps.
   */
  void initMomentum();

  /// Zero total momentum of the system.
  void zeroTotalMomentum();

  /// Return total kinetic energy.
  double kineticEnergy();

  /**
   * Return instantaneous temperature from particle velocities.
   * NOTE: This assumes monotonic, isotropic system for number of degrees of
   * freedom.
   */
  double temperature();

  /// Rescale velocity according to the temperature.
  void rescaleVelocity();

  /// Return status of trial.
  string printStat(const bool header = false);

  /// Update center of mass forces.
  void updateFCOM();

  /// Update half-step velocities.
  void updateVelocityHalfStep();

  /// Perform velocity verlet integration.
  void integrateVelocityVerlet();

  // Return center of mass force for molecule, iMol and dimension, dim.
  double fCOM(const int iMol, const int dim) {
    return fCOM_[space()->dimen()*iMol + dim];
  }

  /// Write restart file.
  void writeRestart(const char* fileName);

  /// Construct from restart file.
  TrialMD(const char* fileName, Pair *pair,
             Criteria *criteria);
  ~TrialMD() {}
  TrialMD* clone(Pair* pair, Criteria* criteria) const {
    TrialMD* t = new TrialMD(*this);
    t->reconstruct(pair, criteria); return t;
  }
  shared_ptr<TrialMD> cloneShrPtr(
    Pair* pair, Criteria* criteria) const {
    return(std::static_pointer_cast<TrialMD, Trial>(
      cloneImpl(pair, criteria)));
  }

 protected:
  vector<double> vel_;    //!< particle velocities
  vector<double> mass_;   //!< particle masses (assumed unity)
  vector<double> fCOM_;   //!< center of mass forces
  string integrator_;     //!< type of integrator

  void attempt1_();

  void defaultConstruction_();

  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl(
    Pair *pair, Criteria *criteria) const {
    shared_ptr<TrialMD> t = make_shared<TrialMD>(*this);
    t->reconstruct(pair, criteria);
    return t;
  }
};

/// Factory method
shared_ptr<TrialMD> makeTrialMD(Pair *pair, Criteria *criteria);

/// Factory method
shared_ptr<TrialMD> makeTrialMD();

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // TRIAL_MD_H_

