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

namespace feasst {

/**
 * Tabular potential for isotropic interacitons.
 */
class PairTabular1D : public Pair {
 public:
  /// Constructor
  explicit PairTabular1D(Space* space, const argtype &args = argtype());

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

 protected:
  vector<vector<double> > rCutInner_;      //!< hard sphere below inner cut-off
  vector<vector<shared_ptr<Table> > > peTable_;  //!< table for potential energy
  string tabFileName_;
  double tol_;      //!< table tolerance

  // See comments of derived class from Pair
  void pairSiteSite_(const int &iSiteType, const int &jSiteType,
    double * energy, double * force, int * neighbor, const double &dx,
    const double &dy, const double &dz);


  // defaults in constructor
  void defaultConstruction_();
};

/// Factory method
shared_ptr<PairTabular1D> makePairTabular1D(Space* space,
  const argtype &args = argtype());

}  // namespace feasst

#endif  // PAIR_TABULAR_1D_H_

