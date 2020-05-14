/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef PAIR_FIELD_SPHERE_H_
#define PAIR_FIELD_SPHERE_H_

#include "./pair_field.h"

namespace feasst {

/**
 * Spherical confinement field potential.
 * Also works as circular confinement in 2D.
 * The distance from the center of the site to the confinement boundary, r,
 * is computed as, \f$r = radius - d\f$, where radius is the radius of the
 * confinement and d is the distance from the center of the site to the center
 * of the confinement.
 */
class PairFieldSphere : virtual public PairField {
 public:
  /// Constructor.
  PairFieldSphere(Space* space);

  /// Initialize the size and orientation of the spherical confinement
  void initSphere(
    /// radius of spherical confinement
    const double radius = 1.,
    /// x-dimension center point of spherical confinement
    const double xSphere = 0.,
    /// y-dimension point
    const double ySphere = 0.,
    /// z-dimension point
    const double zSphere = 0.);

  /// Return the radius
  double radius() const { return radius_; }

  PairFieldSphere(Space* space, const char* fileName);
  ~PairFieldSphere() {}
  virtual PairFieldSphere* clone(Space* space) const {
    PairFieldSphere* p = new PairFieldSphere(*this); p->reconstruct(space); return p;
  }

  void writeRestart(const char* fileName);

 protected:
  double xSphere_, ySphere_, zSphere_, radius_;

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
  virtual void fieldSiteSphere_(
    const int &siteType,  //!< type of site
    double * energy,      //!< energy of interaction
    double * force,       //!< force of interaction
    /// radial distance of center of site to boundary of sphereical confinement
    const double &radial_dist
    ) { ASSERT(0, "not implemented"); }

  void defaultConstruction_();
};

/// Factory method
shared_ptr<PairFieldSphere> makePairFieldSphere(std::shared_ptr<Space> space);

}  // namespace feasst

#endif  // PAIR_FIELD_SPHERE_H_
