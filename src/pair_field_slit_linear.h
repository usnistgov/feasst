/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef PAIR_FIELD_SLIT_LINEAR_H_
#define PAIR_FIELD_SLIT_LINEAR_H_

#include "./pair_field_slit.h"
#include "./pair_field_linear.h"

namespace feasst {

/**
 * Slit confinement field potential with a hard wall and optional square well.
 */
class PairFieldSlitLinear : virtual public PairFieldSlit, virtual public PairFieldLinear {
 public:
  /// Constructor.
  PairFieldSlitLinear(Space* space);

  PairFieldSlitLinear(Space* space, const char* fileName);
  ~PairFieldSlitLinear() {}
  virtual PairFieldSlitLinear* clone(Space* space) const {
    PairFieldSlitLinear* p = new PairFieldSlitLinear(*this); p->reconstruct(space); return p;
  }

  void writeRestart(const char* fileName);

 protected:
  /**
   * Compute the interaction between a site and a field in slit confinement.
   * See base class for documentation of arguments.
   */
  void fieldSiteSlit_(const int &siteType, double * energy, double * force,
    const double &dist) { *energy = linear_energy_(siteType, dist); }

  void defaultConstruction_();
};

/// Factory method
shared_ptr<PairFieldSlitLinear> makePairFieldSlitLinear(std::shared_ptr<Space> space);
shared_ptr<PairFieldSlitLinear> makePairFieldSlitLinear(Space * space);

}  // namespace feasst

#endif  // PAIR_FIELD_SLIT_LINEAR_H_
