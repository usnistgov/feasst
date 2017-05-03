/**
 * \file
 *
 * \brief patchy orientational 'square-well' interactions
 *
 * See the following publications
 *  Fluid–fluid coexistence in colloidal systems with short-ranged strongly
 *  directional attraction
 *  Norbert Kern and Daan Frenkel
 *  J. Chem. Phys. 118, 9882 (2003); http://dx.doi.org/10.1063/1.1569473
 *
 * */
#ifndef PAIR_PATCH_KF_H_
#define PAIR_PATCH_KF_H_

#include "./functions.h"
#include "./pair.h"

namespace feasst {

class PairPatchKF : public Pair {
 public:
  PairPatchKF(Space *space, const double rCut, const double patchAngle);
  ~PairPatchKF();
  virtual PairPatchKF* clone(Space* space) const {
    PairPatchKF* p = new PairPatchKF(*this); p->reconstruct(space); return p;
  }

  int initEnergy();     //!< function to calculate forces, given positions

  /// potential energy of multiple particles
  double multiPartEner(const vector<int> multiPart, const int flag);

  /// potential energy of multiple particles optimized for neighbor list updates
  double multiPartEnerNeigh(const vector<int> multiPart);

  /// stores, restores or updates variables to avoid order recompute of entire
  //  configuration after every change
  void update(const vector<int> mpart, const int flag, const char* uptype);

  double peTot() { return peSR_; }   //!< total potential energy of system

  /// delete particles
  void delPart(const int ipart) {delPartBase(ipart); }
  void delPart(const vector<int> mpart) {delPartBase(mpart); }
  void addPart() {addPartBase(); }                       //!< add one particle

  /// write xyz for visualization
  int printxyz(const char* fileName, const int initFlag,
    const std::string comment="");

  /// sets the mirrorPatch boolean variable
  void mirrorPatch(const int flag) {
    if (flag == 1) { mirrorPatch_ = true; } else { mirrorPatch_ = false; };
  }

  /// read-only access of protected variables
  double peSR() const { return peSR_; }
  double peSRone() const { return peSRone_; }
  double cpa() const { return cpa_; }

 protected:
  const double patchAngle_;   //!< angle of patch (degrees)
  const double cpa_;          //!< cosine of angle of patch
  double peSR_;     //!< total potential energy from lennard-jones interactions
  double peSRone_;  //!< lennard jones potential energy from subset of particles
  double deSR_;               //!< lennard jones potential energy change
  /// mirror patches such that, for each patch, another exists on the other side
  bool mirrorPatch_;
};

}  // namespace feasst

#endif  // PAIR_PATCH_KF_H_

