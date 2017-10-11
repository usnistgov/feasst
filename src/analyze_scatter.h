/**
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of this agreement (see LICENSE.txt) and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#ifndef ANALYZE_SCATTER_H_
#define ANALYZE_SCATTER_H_

#include "./analyze.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Compute scattering, structure factor and radial distribution functions.
 */
class AnalyzeScatter : public Analyze {
 public:
  /// Constructor
  AnalyzeScatter(Pair *pair);

  /// initialize SANS
  void initSANS(const double dgr,  //!< Distance between spatial bins.
    /// Number of macrostates in CriteriaWLTMMC.
    const int nMacros = 1);

  // update analysis every nFreq
  void update() { update(0); }
  void update(const int iMacro);

  // Perform the scattering computation.
  void computeSANS(const int iMacro, const int nMol);
  void computeSANS() { computeSANS(0, space()->nMol()); }

  // print
  void write();
  void write(CriteriaWLTMMC *c);

  // read-only access to protected variables
  vector<vector<vector<vector<long long> > > > histInter() const {
    return histInter2_;
  }
  vector<vector<vector<vector<long long> > > > histIntra() const {
    return histIntra2_;
  }
  vector<double> iq() const { return iq_; }
  vector<vector<double> > iqm() const { return iqm_; }
  vector<double> qwave() const { return qwave_; }
  int qbins() const { return qbins_; }

  /// initialize moments cutoff
  void initMoments(const int nMoments) { nMomentsCut_ = nMoments; }

  AnalyzeScatter(Pair *pair, const char* fileName);
  ~AnalyzeScatter() {}
  AnalyzeScatter* clone(Pair* pair) const {
    AnalyzeScatter* a = new AnalyzeScatter(*this);
    a->reconstruct(pair); return a;
  }
  shared_ptr<AnalyzeScatter> cloneShrPtr(Pair* pair) const {
    return(std::static_pointer_cast<AnalyzeScatter, Analyze>(
      cloneImpl(pair)));
  }
  void writeRestart(const char* fileName) { writeRestart(fileName, -1); };
  void writeRestart(const char* fileName, const int iMacro);

 protected:
  double dgr_;  //!< distance spacing for gr
  int nbins_;   //!< number of bins in gr
  int qbins_;   //!< number of bins for sq
  vector<double> qwave_;  //!< wave vectors
  vector<vector<double> > Pq_;  //!< form factor, Pq[q][itype]
  vector<double> iq_, iqIntra_, iqInter_;   //!< SANS intensity
  vector<vector<double> > iqm_;   //!< SANS intensity for macrostates

  /// intermolecular histogram
  vector<vector<vector<vector<long long> > > > histInter2_;

  /// intramolecular histogram
  vector<vector<vector<vector<long long> > > > histIntra2_;

  vector<long long> countConf_;
  vector<vector<double> > sgr_;   //!< histogram scale factor for g(r)->1

  // Extrapolation quantities
  /// extrapolation histogram
  vector<vector<vector<vector<vector<double> > > > > histMoments_;
  /// maximum number of moments before Taylor series truncation
  ///  e.g., 3 (default) includes moments 0, 1 and 2
  int nMomentsCut_;

  // shared print function
  void printer_(const string fileName, CriteriaWLTMMC *c = NULL,
                const int iMacro = 0);

  /// Return number of particle types from space.
  int nPartTypes_();

  void defaultConstruction_();

  // clone design pattern
  virtual shared_ptr<Analyze> cloneImpl(Pair *pair) const {
    shared_ptr<AnalyzeScatter> a = make_shared<AnalyzeScatter>(*this);
    a->reconstruct(pair); return a;
  }
};

/// Factory method
shared_ptr<AnalyzeScatter> makeAnalyzeScatter(Pair *pair);

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // ANALYZE_SCATTER_H_

