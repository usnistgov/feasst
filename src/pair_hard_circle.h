#ifndef PAIR_HARD_CIRCLE_H_
#define PAIR_HARD_CIRCLE_H_

#include <vector>
#include "./pair.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Hard circles with implicit depletant model (Asakura-Oosawa) based on
 * overlap of excluded area.
 * \f$ U = U_{Hard} + \frac{\Delta A_{ex}}{\pi R_g^2} \phi k_B T\f$
 * The hard particle diameter is set to 1. by default (dCircle_).
 */
class PairHardCircle : public Pair {
 public:
  /// Constructor
  /// @param rCut interaction cut off distance
  PairHardCircle(Space* space, const double rCut);

  /// Initialize the radius of gyration of the depletant, \f$R_g\f$.
  void initRDep(const double rDep = 0.) { rDep_ = rDep; }

  /// Initialize interactions.
  virtual void initEnergy();

  /// Potential energy of multiple particles
  double multiPartEner(const vector<int> multiPart, const int flag);
  void multiPartEnerAtomCutInner(const double &r2, const int &itype,
    const int &jtype);

/**
 * Potential energy and forces of all particles
 *  if flag == 0, dummy calculation
 *  if flag == 1, all config calculation
 */
  double allPartEnerForce(const int flag);

  /// potential energy and forces of all particles
  double allPartEnerForceNoCell();

  /// Stores, restores or updates variables to avoid order recompute of entire
  ///  configuration after every change.
  void update(const vector<int> mpart, const int flag, const char* uptype);

  /// read-only access of protected variables
  double peSR() const { return peSR_; }
  double peSRone() const { return peSRone_; }

  /// Write restart file.
  virtual void writeRestart(const char* fileName);

  /// Construct from restart file.
  PairHardCircle(Space* space, const char* fileName);
  virtual ~PairHardCircle() {}
  virtual PairHardCircle* clone(Space* space) const {
    PairHardCircle* p = new PairHardCircle(*this);
    p->reconstruct(space);
    return p;
  }

 protected:
  double peSR_;
  double deSR_;             //!< potential energy change
  double dCircle_;               //!< diameter of hard circle
  double rDep_;               //!< radius of depletant

  // defaults in constructor
  void defaultConstruction_();
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // PAIR_HARD_CIRCLE_H_

