/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef TRIAL_CLUSTER_H_
#define TRIAL_CLUSTER_H_

#include <memory>
#include <string>
#include "./trial.h"

namespace feasst {

/**
 * Attempt a rigid cluster rotation or translation of all clusters in system.
 */
class TrialCluster : public Trial {
 public:
  /// Constructor.
  /// @param transType For rigid translation, use "clustertrans".
  ///   For rigid rotations, use "clusterrotate".
  TrialCluster(Pair *pair, Criteria *criteria,
               const char* transType);

  /// Cut off for distance-based cluster criteria.
  /// If the particles centers within clusterCut, they are in a cluster.
  double clusterCut;

  /// This constructor is not often used, but its purpose is to initialize trial
  /// for interface before using reconstruct to set object pointers.
  explicit TrialCluster(const char* transType);

  /// Write restart file.
  void writeRestart(const char* fileName);

  /// Construct from restart file.
  TrialCluster(const char* fileName, Pair *pair,
               Criteria *criteria);

  ~TrialCluster() {}
  TrialCluster* clone(Pair* pair, Criteria* criteria) const {
    TrialCluster* t = new TrialCluster(*this);
    t->reconstruct(pair, criteria); return t; }
  shared_ptr<TrialCluster> cloneShrPtr
    (Pair* pair, Criteria* criteria) const {
      return(std::static_pointer_cast<TrialCluster, Trial>
      (cloneImpl(pair, criteria))); }

  /// tune max parameters for translation and rotation moves
  void tuneParameters();

  double targAcceptPer;      //!< target acceptance percentage

  // functions for read-only access of private data-members
  string transType() const { return transType_; }

 protected:
  string transType_;  //!< type of transformation

  void attempt1_();

  void defaultConstruction_();

  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl
    (Pair *pair, Criteria *criteria) const {
      shared_ptr<TrialCluster> t = make_shared<TrialCluster>(*this);
      t->reconstruct(pair, criteria);
      return t;
  }
};

shared_ptr<TrialCluster> makeTrialCluster(Pair *pair,
  Criteria *criteria, const char* transType);

shared_ptr<TrialCluster> makeTrialCluster(const char* transType);

class MC;

/// Add a "TrialCluster" object to Monte Carlo object, mc
void clusterTrial(MC *mc, const char* transType, const double clusterCut = -1,
                  const double maxMoveParam = -1);

/// Add a "TrialCluster" object to Monte Carlo object, mc
void clusterTrial(shared_ptr<MC> mc, const char* type,
                  const double clusterCut = -1, const double maxMoveParam = -1);

}  // namespace feasst

#endif  // TRIAL_CLUSTER_H_

