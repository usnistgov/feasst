/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef PAIR_FIELD_LINEAR_H_
#define PAIR_FIELD_LINEAR_H_

#include "./pair.h"

namespace feasst {

/**
 * Interaction with a piece of the confinement is given by a hard potential
 * and an optional attraction within a distance of the hard surface,
 *
 * \f$ U = \begin{cases}
 *         \infty               & r < \sigma/2
 *         \epsilon (r\rCut-1)  & \sigma/2 \leq r \leq rCut
 *         0                    & otherwise
 *         \end{cases}\f$
 *
 * where r is the distance from the center of the site to the confinement, and
 * rCut is may be initialied by Pair::initRCut.
 *
 * Note that a positive epsilon leads to an attraction.
 *
 * For hard walls only, use Pair::halfSig2rCut to remove attractions by setting rCut
 * to \f$\sigma/2\f$.
 */
class PairFieldLinear : virtual public Pair {
 public:
  /// Constructor
  PairFieldLinear(Space * space) : Pair(space) { defaultConstruction_(); }

  PairFieldLinear(Space * space, const char* fileName);
  void writeRestart(const char* fileName);

 protected:
  void defaultConstruction_();

  /// Return the energy of interaction for a distance, dist, from confinement
  double linear_energy_(const int &siteType, const double &dist);
};

}  // namespace feasst

#endif  // PAIR_FIELD_LINEAR_H_
