/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef PAIR_FIELD_LJ_H_
#define PAIR_FIELD_LJ_H_

#include "./pair.h"

namespace feasst {

/**
 * Interaction with an infinitesimal piece of the confinement is given by a
 * power law,
 *
 * \f$dU \approx \epsilon \left( \frac{\sigma}{r + \Delta} \right)^\alpha\f$
 *
 * where r is the distance from the center of the site to the confinement.
 * Note that the scaling for the total potential energy depends upon
 * integration over the confinement geometry.
 * \f$\epsilon\f$ and \f$\sigma\f$ are set by the data file or Pair::initEps
 * and Pair::initSig.
 */
class PairFieldLJ : virtual public Pair {
 public:
  /// Constructor
  PairFieldLJ(Space * space) : Pair(space) { defaultConstruction_(); }

  /// Initialize the \f$\alpha\f$ exponent
  void initAlpha(const double alpha = 6.) { alpha_ = alpha; }

  /// Initialize the \f$\Delta\f$ shift factor
  void initDelta(const double delta = 0.) { delta_ = delta; }

  PairFieldLJ(Space * space, const char* fileName);
  void writeRestart(const char* fileName);

 protected:
  double alpha_, delta_;
  void defaultConstruction_();
};

}  // namespace feasst

#endif  // PAIR_FIELD_LJ_H_
