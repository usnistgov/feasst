/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef PAIR_PATCH_KF_H_
#define PAIR_PATCH_KF_H_

#include <string>
#include <vector>
#include "./functions.h"
#include "./pair.h"

namespace feasst {

/**
 * Patchy orientational dependent square-well interactions.
 * See the following publications:
 *  Fluid–fluid coexistence in colloidal systems with short-ranged strongly
 *  directional attraction
 *  Norbert Kern and Daan Frenkel
 *  J. Chem. Phys. 118, 9882 (2003); http://dx.doi.org/10.1063/1.1569473
 *
 * As currently implemented, this class assumes only one central particle site
 * with one or more patches, and not multiple particles, each with their own
 * set of patches, bonded together.
 */
class PairPatchKF : public Pair {
 public:
  /// Constructor
  PairPatchKF(Space *space,
   /**
    * allowed string key pairs (e.g., dictionary):
    *
    * patchAngle : solid half-angle (in degrees) of patch (see reference)
    */
   const argtype &args = argtype());

  /// Mirrors the patch on the other side of the particle if flag == 1.
  void mirrorPatch(const int flag) {
    if (flag == 1) { mirrorPatch_ = true; } else { mirrorPatch_ = false; };
  }

  /// potential energy of multiple particles
  double multiPartEner(const vector<int> multiPart, const int flag);

  /// potential energy of multiple particles optimized for neighbor list updates
  virtual double multiPartEnerNeigh(const vector<int> multiPart);

  /// stores, restores or updates variables to avoid order recompute of entire
  // configuration after every change
  void update(const vector<int> mpart, const int flag, const char* uptype);

  /// Write xyz for visualization.
  /// @param initFlag open if flag is 1, append if flag is 0.
  virtual int printXYZ(const char* fileName, const int initFlag,
    const std::string comment = "");

  /// Update clusters of entire system.
  /// @param tol unused parameter.
  void updateClusters(const double tol);

  /// potential energy and forces of all particles
  virtual double allPartEnerForce(const int flag);

  /// Compute interactions between two molecules.
  virtual void allPartEnerForceMolCutInner(const double r2,
    const int iMol, const int jMol, const double dx,
    const double dy, const double dz);

  /// Return the cosine of the patch angle.
  double cpa() const { return cpa_; }

  PairPatchKF(Space* space, const char* fileName) : Pair(space) {
    ASSERT(0, "no restart implemented"); }
  ~PairPatchKF() {};
  virtual PairPatchKF* clone(Space* space) const {
    PairPatchKF* p = new PairPatchKF(*this); p->reconstruct(space); return p;
  }

 protected:
  double patchAngle_;   //!< angle of patch (degrees)
  double cpa_;          //!< cosine of angle of patch
  double deSR_;               //!< lennard jones potential energy change
  /// mirror patches such that, for each patch, another exists on the other side
  bool mirrorPatch_;
};

/// Factory method
shared_ptr<PairPatchKF> makePairPatchKF(Space *space, const argtype &args);

}  // namespace feasst

#endif  // PAIR_PATCH_KF_H_

