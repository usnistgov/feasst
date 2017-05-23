/**
 * \file
 *
 * \brief
 *
 */

#ifndef ANALYZE_SCATTER_H_
#define ANALYZE_SCATTER_H_

#include "./analyze.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

class AnalyzeScatter : public Analyze {
 public:
  AnalyzeScatter(Space *space, Pair *pair);
  AnalyzeScatter(Space *space, Pair *pair, const char* fileName);
  ~AnalyzeScatter() {}
  AnalyzeScatter* clone(Space* space, Pair* pair) const {
    AnalyzeScatter* a = new AnalyzeScatter(*this);
    a->reconstruct(space, pair); return a;
  }
  shared_ptr<AnalyzeScatter> cloneShrPtr(Space* space, Pair* pair) const {
    return(std::static_pointer_cast<AnalyzeScatter, Analyze>(
      cloneImpl(space, pair)));
  }
  void defaultConstruction();
  void writeRestart(const char* fileName, const int iMacro = -1);

  /// initialize SANS
  void initSANS(const double dgr, const int nMacros);
  void initSANS(const double dgr) { initSANS(dgr, 1); }

  /// update analysis every nFreq
  void update() { update(0); }
  void update(const int iMacro);

  /// compute
  void computeSANS(const int iMacro, const int nMol);
  void computeSANS() { computeSANS(0, space_->nMol()); }

  /// print
  void write();
  void write(CriteriaWLTMMC *c);

  /// determine number of particle types
  int nPartTypes();

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

  // shared print function
  void printer_(const string fileName, CriteriaWLTMMC *c = NULL,
                const int iMacro = 0);

  // clone design pattern
  virtual shared_ptr<Analyze> cloneImpl(Space* space, Pair *pair) const {
    shared_ptr<AnalyzeScatter> a = make_shared<AnalyzeScatter>(*this);
    a->reconstruct(space, pair); return a;
  }
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // ANALYZE_SCATTER_H_

