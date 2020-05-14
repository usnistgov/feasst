/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef PAIR_FIELD_SLIT_H_
#define PAIR_FIELD_SLIT_H_

#include "./pair_field.h"

namespace feasst {

/**
 * Fields in slit confinement depend upon the identity of the site and r, the
   nearest distance from the center of the site to the upper and lower
   boundaries.

   The distance, r, is computed for both the upper and lower boundary.

   \f$ r = \begin{cases}
           upper - x \\
           x - lower
           \end{cases}\f$

   where x is the position of the center of the site in the dimension of the
   slit confinement.
 * A negative r indicates the particle is outside of the slit confinement.
 */
class PairFieldSlit : virtual public PairField {
 public:
  /// Constructor.
  PairFieldSlit(Space* space);

  /// Initialize the size and orientation of the slit confinement
  void initSlit(
    /// dimension index which is perpendicular to the interface of the slit confinement
    const int dimension = -1,
    /// upper bound of the slit confinement
    const double upper = 0,
    /// lower bound of the slit confinement
    const double lower = 0);

  PairFieldSlit(Space* space, const char* fileName);
  ~PairFieldSlit() {}
  virtual PairFieldSlit* clone(Space* space) const {
    PairFieldSlit* p = new PairFieldSlit(*this); p->reconstruct(space); return p;
  }

  void writeRestart(const char* fileName);

 protected:
  double upper_, lower_;
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
   * Compute the interaction between a site and a field in a slit confinement.
   */
  virtual void fieldSiteSlit_(
    const int &siteType,  //!< type of site
    double * energy,      //!< energy of interaction
    double * force,       //!< force of interaction
    const double &dist    //!< distance from center of site to slit
    ) { ASSERT(0, "not implemented"); }

  void defaultConstruction_();
};

/// Factory method
shared_ptr<PairFieldSlit> makePairFieldSlit(std::shared_ptr<Space> space);

}  // namespace feasst

#endif  // PAIR_FIELD_SLIT_H_
