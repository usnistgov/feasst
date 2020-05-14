/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef PAIR_FIELD_SLIT_LJ_H_
#define PAIR_FIELD_SLIT_LJ_H_

#include "./pair_field_slit.h"
#include "./pair_field_lj.h"

namespace feasst {

/**
 * Slit field potential with a power law and optional position shift,
 *
 * \f$U = \begin{cases}
 *        \infty & r \leq 0 \\
 *        \epsilon \left( \sigma / r \right)^\alpha & otherwise
 *        \end{cases}\f$
 *
 * where r is the distance from the bounadries to the center of the site.
 * Note that r is defined precisely in PairFieldSlit.
 */
class PairFieldSlitLJ : virtual public PairFieldSlit,
  virtual public PairFieldLJ {
 public:
  /// Constructor.
  PairFieldSlitLJ(Space* space);
  PairFieldSlitLJ(Space* space, const char* fileName);
  ~PairFieldSlitLJ() {}
  virtual PairFieldSlitLJ* clone(Space* space) const {
    PairFieldSlitLJ* p = new PairFieldSlitLJ(*this); p->reconstruct(space); return p;
  }

  void writeRestart(const char* fileName);

 protected:
  /**
   * Compute the interaction between a site and a field in slit confinement.
   * See base class for documentation of arguments.
   */
  void fieldSiteSlit_(const int &siteType, double * energy, double * force,
    const double &dist);

  void defaultConstruction_();
};

/// Factory method
shared_ptr<PairFieldSlitLJ> makePairFieldSlitLJ(std::shared_ptr<Space> space);

}  // namespace feasst

#endif  // PAIR_FIELD_SLIT_LJ_H_
