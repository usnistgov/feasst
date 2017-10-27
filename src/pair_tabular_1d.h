/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef PAIR_TABULAR_1D_H_
#define PAIR_TABULAR_1D_H_

#include <memory>
#include <string>
#include <vector>
#include "./pair.h"
#include "./table.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Tabular potential for isotropic interacitons.
 */
class PairTabular1D : public Pair {
 public:
  /// Constructor
  explicit PairTabular1D(Space* space);

  /// Read table from file.
  void readTable(const char* fileName);

  /// Set the interpolator.
  void setInterpolator(const char* name);

  /// Read-only access to tables.
  vector<vector<shared_ptr<Table> > > peTable() const { return peTable_; }

  /// Write restart file.
  virtual void writeRestart(const char* fileName);

  /// Construct from restart file.
  PairTabular1D(Space* space, const char* fileName);
  virtual ~PairTabular1D() {}
  virtual PairTabular1D* clone(Space* space) const;

  /// potential energy and forces of all particles
  void multiPartEnerAtomCutInner(const double &r2, const int &itype,
    const int &jtype);

 protected:
  vector<vector<double> > rCutInner_;      //!< hard sphere below inner cut-off
  vector<vector<shared_ptr<Table> > > peTable_;  //!< table for potential energy
  string tabFileName_;
  double tol_;      //!< table tolerance

  // defaults in constructor
  void defaultConstruction_();
};

/// Factory method
shared_ptr<PairTabular1D> makePairTabular1D(Space* space);

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // PAIR_TABULAR_1D_H_

