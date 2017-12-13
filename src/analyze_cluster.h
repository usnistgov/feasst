/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef ANALYZECLUSTER_H_
#define ANALYZECLUSTER_H_

#include "./analyze.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Compute cluster properties
 */
class AnalyzeCluster : public Analyze {
 public:
  /// Constructor
  AnalyzeCluster(Pair *pair, const argtype &args = argtype());

  /// Initialize distance-based cluster cutoff.
  void initClusterCut(const double clusterCut) { clusterCut_ = clusterCut; }

  /// Initialize size of block average.
  void initBlockAverage(const int nBlock) { nClusterAccVec_.setBlock(nBlock); }

  // update analysis every nFreq
  void update() { update(0); }
  void update(const int iMacro);

  // print
  void write(CriteriaWLTMMC *c);

  double zbin;  // histogram bin size

  /**
   * Compute percolation probability.
   * A percolating cluster is one in which a particle in the cluster is
   * connected to its own periodic image.
   */
  void initPercolation(
    const int percFlag = 0
    /**< if "0", no computation.
         if "1", use expanding box (slow).
         if "2", use contact map. */
    ) { percFlag_ = percFlag; }

  /// Return the number of clusters
  AccumulatorVec nClusterAccVec() const { return nClusterAccVec_; }

  /// Return the coordination number
  AccumulatorVec coordNumAccVec() const { return coordNumAccVec_; }

  /// Return the percolation
  AccumulatorVec percolation() const { return percolation_; }

  /// Write restart file.
  void writeRestart(const char* fileName);

  /// Construct from restart file.
  AnalyzeCluster(Pair *pair, const char* fileName);

  ~AnalyzeCluster() {}
  AnalyzeCluster* clone(Pair* pair) const {
    AnalyzeCluster* a = new AnalyzeCluster(*this);
    a->reconstruct(pair); return a;
  }
  shared_ptr<AnalyzeCluster> cloneShrPtr(Pair* pair) const {
    return(std::static_pointer_cast<AnalyzeCluster, Analyze>(cloneImpl(
      pair)));
  }

  // Nematic matrix formed by average outer product of directors
  // https://doi.org/10.1103/PhysRevLett.98.108303
  // Nematic, S = max(eigen(1/2 <3nematic-I>))
  vector<vector<AccumulatorVec> > nematic_;

  /// Instantaneous (single config) coordination number for a given molecule
  vector<int> iMol2coord() const { return iMol2coord_; }

 protected:
  double clusterCut_;                 // cluster cut-off definition
  AccumulatorVec nClusterAccVec_;     // number of clusters in system
  AccumulatorVec coordNumAccVec_;     // average coordination number
  AccumulatorVec largestClusAccVec_;  // average largest cluster
  AccumulatorVec percolation_;        // 0 or 1 if system spanning
  AccumulatorVec pcos_;

  int percFlag_;   //!< type of percolation computation

  // contact orientation
  vector<Histogram> zOrient_;     // histogram of z-axis orientation

  vector<Histogram> nClusSize_;   // histogram of cluster size

  vector<int> iMol2coord_;

  void defaultConstruction_();

  // clone design pattern
  virtual shared_ptr<Analyze> cloneImpl(Pair *pair) const {
    shared_ptr<AnalyzeCluster> a = make_shared<AnalyzeCluster>(*this);
    a->reconstruct(pair); return a;
  }
};

/// Factory method
shared_ptr<AnalyzeCluster> makeAnalyzeCluster(Pair *pair,
  const argtype &args = argtype());

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // ANALYZECLUSTER_H_

