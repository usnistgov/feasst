/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef PAIR_FIELD_CYLINDER_SW_H_
#define PAIR_FIELD_CYLINDER_SW_H_

#include "./pair_field_cylinder.h"
#include "./pair_field_sw.h"

namespace feasst {

/**
 * Cylinder confinement field potential with a hard wall and optional square well.
 */
class PairFieldCylinderSW : virtual public PairFieldCylinder,
  virtual public PairFieldSW {
 public:
  /// Constructor.
  PairFieldCylinderSW(Space* space);

  PairFieldCylinderSW(Space* space, const char* fileName);
  ~PairFieldCylinderSW() {}
  virtual PairFieldCylinderSW* clone(Space* space) const {
    PairFieldCylinderSW* p = new PairFieldCylinderSW(*this); p->reconstruct(space); return p;
  }

  void writeRestart(const char* fileName);

 protected:
  /**
   * Compute the interaction between a site and a field in spherical
   * confinement.
   * See base class for documentation of arguments.
   */
  void fieldSiteCylinder_(const int &siteType, double * energy, double * force,
    const double &radial_dist) { *energy += sw_energy_(siteType, radial_dist); }

  void defaultConstruction_();
};

/// Factory method
shared_ptr<PairFieldCylinderSW> makePairFieldCylinderSW(std::shared_ptr<Space> space);

}  // namespace feasst

#endif  // PAIR_FIELD_CYLINDER_SW_H_
