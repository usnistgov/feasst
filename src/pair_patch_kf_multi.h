/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef PAIR_PATCH_KF_MULTI_H_
#define PAIR_PATCH_KF_MULTI_H_

#include <string>
#include <vector>
#include "./functions.h"
#include "./pair_patch_kf.h"
#include "./bond.h"

namespace feasst {

/**
 * Patchy orientational dependent square-well interactions.
 *
 * See the following publications:
 *  Fluid–fluid coexistence in colloidal systems with short-ranged strongly
 *  directional attraction
 *  Norbert Kern and Daan Frenkel
 *  J. Chem. Phys. 118, 9882 (2003); http://dx.doi.org/10.1063/1.1569473
 *
 * A patch interaction is defined in a LAMMPS-like data file by two sites,
 * the patch center and the patch orienter.
 * The patch center is the center of a particle with excluded volume
 * (e.g., the hard sphere in a Kern-Frenkel patch).
 * The unit vector which describes the orientation of the patch is given by
 * the difference between the patch orienter and the patch center.
 * The patch center must have a non-zero sigma parameter, the second pair
 * coefficient.
 * The patch orienter must have a zero sigma parameter, but the epsilon
 * is the strength of the patch interaction.
 * The patch orienter must also be bonded to the patch center.
 */
class PairPatchKFMulti : public PairPatchKF {
 public:
  /// Constructor
  PairPatchKFMulti(Space *space,
   /**
    * allowed string key pairs (e.g., dictionary):
    *
    * patchAngle : Solid half-angle (in degrees) of patch (see reference).
    *              Note the patch angle is from the center to edge of cone.
    *              All particle types use this patch angle.
    */
   const argtype &args = argtype());

  /// Initialize potential cutoff (square-well distance for patches)
  /// based on the given value of rCut for all beads.
  /// Also sets patch angle the same for all patches.
  void initIJ();

  // HWH: Depreciate: use initIJ()
  void initCuts() { initIJ(); }

  // overload virtual to test bond parameters
  void initData(const string fileName);

  /// Set the patch angle.
  void initPatchAngle(const double angle);

  /// Initialize the patch angle
  void initPatchAngleInDegrees(const double angle, const int itype);

  /// Magnitude of pair interaction must be larger than this to be a clsuter
  double clusterTolerance = 0.;

  /// potential energy of multiple particles optimized for neighbor list updates
  double multiPartEnerNeigh(const vector<int> multiPart);

  /// Write xyz for visualization.
  /// @param initFlag open if flag is 1, append if flag is 0.
  int printXYZ(const char* fileName, const int initFlag,
    const std::string comment = "");

  /// potential energy and forces of all particles
  double allPartEnerForce(const int flag);

  // ANN: customize patch function here
  /// Patch type
  /// if 0, kern frenkel
  /// if 1, https://doi.org/10.1063/1.5021347
  int patchType = 0;

  PairPatchKFMulti(Space* space, const char* fileName) : PairPatchKF(space) {
    ASSERT(0, "no restart implemented"); }
  ~PairPatchKFMulti() {};
  virtual PairPatchKFMulti* clone(Space* space) const {
    PairPatchKFMulti* p = new PairPatchKFMulti(*this); p->reconstruct(space); return p;
  }

  // overloaded virtual function
  void reconstruct(Space* space);

  /// Set the order parameter value.
  virtual void setOrder(const double order);

 protected:
  vector<double> cpaSqi_; // cosine(patchAngle)^2 for a given type

  shared_ptr<Bond> bond_;
};

/// Factory method
shared_ptr<PairPatchKFMulti> makePairPatchKFMulti(Space *space, const argtype &args);

}  // namespace feasst

#endif  // PAIR_PATCH_KF_MULTI_H_

