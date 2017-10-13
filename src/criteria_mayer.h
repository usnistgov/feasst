/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#ifndef CRITERIA_MAYER_H_
#define CRITERIA_MAYER_H_

#include "./criteria.h"
#include "./pair.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Mayer Sampling: Calculation of Cluster Integrals using Free-Energy
 * Perturbation Methods
 * Jayant K. Singh and David A. Kofke, Phys. Rev. Lett., 92, 220601, 2004.
 * https://doi.org/10.1103/PhysRevLett.92.220601
 */
class CriteriaMayer : public Criteria {
 public:
  /// Constructor
  explicit CriteriaMayer(const double beta);

  /// Construct by checkpoint file.
  explicit CriteriaMayer(const char* fileName);

  ~CriteriaMayer() {}
  CriteriaMayer* clone() const;
  shared_ptr<CriteriaMayer> cloneShrPtr() const;

  /// Return whether to accept (1) or reject (0).
  int accept(const double lnpMet, const double peNew, const char* moveType,
    const int reject);

  /// Store reference and full potential energies of the current configuration.
  void store(Pair* pair);

  /// Initialize reference potential.
  void initPairRef(Pair *pair) { pairRef_ = pair; }

  /// Mayer sampling ensemble averages of full potential.
  Accumulator mayer;

  /// Mayer sampling ensemble averages of reference potential.
  Accumulator mayerRef;

  /// Return ratio of virial coefficients of the full potential by reference.
  double b2ratio() const { return mayer.average()/mayerRef.average(); };

 protected:
  Pair* pairRef_;

  /// configuration of old configuration, recorded by "store"
  double peOld_, peRefOld_;

  /// defaults in constructor
  void defaultConstruction_();

  // clone design pattern
  virtual shared_ptr<Criteria> cloneImpl_() const;
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // CRITERIA_MAYER_H_

