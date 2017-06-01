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

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

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

  /// write xyz for visualization
  int printxyz(const char* fileName, const int initFlag,
    const std::string comment="");

  /// sets the mirrorPatch boolean variable
  void mirrorPatch(const int flag) {
    if (flag == 1) { mirrorPatch_ = true; } else { mirrorPatch_ = false; };
  }

  /// update clusters of entire system
  void updateClusters(const double tol  //!< unused parameter
    );

  /// potential energy and forces of all particles
  double allPartEnerForce(const int flag);
  
  /// read-only access of protected variables
  double cpa() const { return cpa_; }

 protected:
  const double patchAngle_;   //!< angle of patch (degrees)
  const double cpa_;          //!< cosine of angle of patch
  double deSR_;               //!< lennard jones potential energy change
  /// mirror patches such that, for each patch, another exists on the other side
  bool mirrorPatch_;
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // PAIR_PATCH_KF_H_

