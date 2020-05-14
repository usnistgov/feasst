/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef PAIR_FIELD_CYLINDER_H_
#define PAIR_FIELD_CYLINDER_H_

#include "./pair_field.h"

namespace feasst {

/**
 * Cylinderical confinement field potential.
 * The distance from the center of the site to the confinement boundary, r,
 * is computed as, \f$r = radius - d\f$, where radius is the radius of the
 * confinement and d is the nearest distance from the center of the site to the
 * axis of the confinement.
 */
class PairFieldCylinder : virtual public PairField {
 public:
  /// Constructor.
  PairFieldCylinder(Space* space);

  /// Initialize the size and orientation of the spherical confinement
  void initCylinder(
    /// Cartesian dimension that the cylindrical axis lies along
    const int dimen = 0,
    /// radius of cylindrical confinement
    const double radius = 1.,
    /// x-dimension of any point along axis of cylindrical confinement
    const double xCylinder = 0.,
    /// y-dimension point
    const double yCylinder = 0.,
    /// z-dimension point
    const double zCylinder = 0.);

  /// Return the radius
  double radius() const { return radius_; }

  PairFieldCylinder(Space* space, const char* fileName);
  ~PairFieldCylinder() {}
  virtual PairFieldCylinder* clone(Space* space) const {
    PairFieldCylinder* p = new PairFieldCylinder(*this); p->reconstruct(space); return p;
  }

  void writeRestart(const char* fileName);

 protected:
  double xCylinder_, yCylinder_, zCylinder_, radius_;
  int confine_dimen_;

  /**
   * Compute the interaction between a site and a field.
   */
  void fieldSite_(
    const int &siteType,  //!< type of site
    double * energy,      //!< energy of interaction
    double * force,       //!< force of interaction
    const double &x,      //!< x-dimension absolute position
    const double &y,      //!< y-dimension absolute position
    const double &z);     //!< z-dimension absolute position

  /**
   * Compute the interaction between a site and a field in spherical confinement.
   */
  virtual void fieldSiteCylinder_(
    const int &siteType,  //!< type of site
    double * energy,      //!< energy of interaction
    double * force,       //!< force of interaction
    /// radial distance of center of site to cylindrical confinement boundary
    const double &radial_dist
    ) { ASSERT(0, "not implemented"); }

  void defaultConstruction_();
};

/// Factory method
shared_ptr<PairFieldCylinder> makePairFieldCylinder(std::shared_ptr<Space> space);

}  // namespace feasst

#endif  // PAIR_FIELD_CYLINDER_H_
