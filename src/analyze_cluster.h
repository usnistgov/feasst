/**
 * \file
 *
 * \brief
 *
 */

#ifndef ANALYZECLUSTER_H_
#define ANALYZECLUSTER_H_

#include "./analyze.h"

namespace feasst {

class AnalyzeCluster : public Analyze {
 public:
  AnalyzeCluster(Space *space, Pair *pair);
  AnalyzeCluster(Space *space, Pair *pair, const char* fileName);
  ~AnalyzeCluster() {}
  AnalyzeCluster* clone(Space* space, Pair* pair) const {
    AnalyzeCluster* a = new AnalyzeCluster(*this);
    a->reconstruct(space, pair); return a;
  }
  shared_ptr<AnalyzeCluster> cloneShrPtr(Space* space, Pair* pair) const {
    return(std::static_pointer_cast<AnalyzeCluster, Analyze>(cloneImpl(
      space, pair)));
  }
  void defaultConstruction();
  void writeRestart(const char* fileName);

  /// initialize cluster cutoff
  void initClusterCut(const double clusterCut) { clusterCut_ = clusterCut; }

  /// initialize size of block average
  void initBlockAverage(const int nBlock) { nClusterAccVec_.setBlock(nBlock); }

  /// update analysis every nFreq
  void update() { update(0); }
  void update(const int iMacro);

  /// print
  void print(CriteriaWLTMMC *c);

  double zbin;                   // histogram bin size

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

  /// read-only access
  AccumulatorVec nClusterAccVec() const { return nClusterAccVec_; }
  AccumulatorVec coordNumAccVec() const { return coordNumAccVec_; }
  AccumulatorVec percolation() const { return percolation_; }

 protected:
  double clusterCut_;                 // cluster cut-off definition
  AccumulatorVec nClusterAccVec_;     // number of clusters in system
  AccumulatorVec coordNumAccVec_;     // average coordination number
  AccumulatorVec largestClusAccVec_;  // average largest cluster
  AccumulatorVec percolation_;        // 0 or 1 if system spanning
  
  int percFlag_;   //!< type of percolation computation

  // contact orientation
  vector<Histogram> zOrient_;     // histogram of z-axis orientation

  vector<Histogram> nClusSize_;   // histogram of cluster size

  // clone design pattern
  virtual shared_ptr<Analyze> cloneImpl(Space* space, Pair *pair) const {
    shared_ptr<AnalyzeCluster> a = make_shared<AnalyzeCluster>(*this);
    a->reconstruct(space, pair); return a;
  }
};

}  // namespace feasst

#endif  // ANALYZECLUSTER_H_

