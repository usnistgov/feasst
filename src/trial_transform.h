/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef TRIAL_TRANSFORM_H_
#define TRIAL_TRANSFORM_H_

#include <memory>
#include <string>
#include "./trial.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Attempt a transformation of particle position(s) or simulation box.
 */
class TrialTransform : public Trial {
 public:
  /**
   * Constructor
   * @param transType
   *  For rigid single-particle translations, "translate".
   *  For rigid single-particle rotations, "rotate".
   *  For "smart Monte Carlo, force biased", use "smctrans".
   *    Note that "smctrans" is for translational moves only, and not rotation,
   *    and these require center-of-mass forces.
   *    See http://dx.doi.org/10.1063/1.436415
   *  For triclinic cells, "xyztilt", "xztilt" and "yztilt".
   *  For volume change, "vol"
   *  For x-dimension box length change, "lxmod". Similarly "lymod", "lzmod".
   */
  TrialTransform(Pair *pair, Criteria *criteria,
                 const char* transType);

  /// Constructor
  TrialTransform(Pair *pair, Criteria *criteria,
    const argtype &args = argtype());

  /// This constructor is not often used, but its purpose is to initialize trial
  /// for interface before using reconstruct to set object pointers.
  explicit TrialTransform(const char* transType);

  // tune parameters (e.g., based on acceptance)
  void tuneParameters();
  double targAcceptPer;      //!< target acceptance percentage

  /// accumulator for order parameter of trial
  Accumulator paramAccumulator;

  /// Select molecule/particle type for transformations.
  /// Currently only implemented for translate/rotate.
  void selectType(const char* molType) { molType_ = molType; }

  /// Write restart file.
  void writeRestart(const char* fileName);

  /// Construct from restart file.
  TrialTransform(const char* fileName, Pair *pair,
                 Criteria *criteria);
  ~TrialTransform() {}
  TrialTransform* clone(Pair* pair, Criteria* criteria) const {
    TrialTransform* t = new TrialTransform(*this);
    t->reconstruct(pair, criteria); return t;
  }
  shared_ptr<TrialTransform> cloneShrPtr(
    Pair* pair, Criteria* criteria) const {
    return(std::static_pointer_cast<TrialTransform, Trial>(
      cloneImpl(pair, criteria)));
  }

  /// Return status of trial.
  string printStat(const bool header = false);

  /// Return transType.
  string transType() const { return transType_; }

 protected:
  string transType_;  //!< type of transformation
  string molType_;    //!< type of molecule to transform

  void attempt1_();

  /// Attempt to scale the domain by a factor.
  void scaleAttempt_(const double factor);

  void defaultConstruction_();

  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl(
    Pair *pair, Criteria *criteria) const {
    shared_ptr<TrialTransform> t = make_shared<TrialTransform>(*this);
    t->reconstruct(pair, criteria);
    return t;
  }
};

/// Factory method
shared_ptr<TrialTransform> makeTrialTransform(Pair *pair,
  Criteria *criteria, const char* transType);

/// Factory method
shared_ptr<TrialTransform> makeTrialTransform(Pair *pair, Criteria *criteria,
  const argtype &args = argtype());

/// Factory method
shared_ptr<TrialTransform> makeTrialTransform(const char* transType);

class MC;

/// Add a "TrialTransform" object to the Monte Carlo object, mc.
void transformTrial(MC *mc, const char* transType, double maxMoveParam = -1);

/// Add a "TrialTransform" object to the Monte Carlo object, mc.
void addTrialTransform(MC *mc, const argtype &args = argtype());

// Renaming of above:
void transformTrial(MC *mc, const argtype &args = argtype());

/// Add a "TrialTransform" object to the Monte Carlo object, mc.
void transformTrial(shared_ptr<MC> mc, const char* transType,
                    double maxMoveParam = -1);

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // TRIAL_TRANSFORM_H_

