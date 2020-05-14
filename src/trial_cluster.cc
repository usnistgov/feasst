/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./trial_cluster.h"
#include "mc.h"

namespace feasst {

TrialCluster::TrialCluster(
  const char* transType)
  : Trial(),
    transType_(transType) {
  defaultConstruction_();
}
TrialCluster::TrialCluster(
  Pair *pair,
  Criteria *criteria,
  const char* transType)
  : Trial(pair, criteria),
    transType_(transType) {
  defaultConstruction_();
}
TrialCluster::TrialCluster(const char* fileName,
  Pair *pair,
  Criteria *criteria)
  : Trial(pair, criteria, fileName) {
  transType_ = fstos("transType", fileName);
  defaultConstruction_();
  maxMoveParam = fstod("maxMoveParam", fileName);
  targAcceptPer = fstod("targAcceptPer", fileName);
  clusterCut = fstod("clusterCut", fileName);
}

void TrialCluster::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  file << "# transType " << transType_ << endl;
  file << "# maxMoveParam " << maxMoveParam << endl;
  file << "# targAcceptPer " << targAcceptPer << endl;
  file << "# clusterCut " << clusterCut << endl;
}

void TrialCluster::defaultConstruction_() {
  className_.assign("TrialCluster");
  trialType_.assign("move");
  verbose_ = 0;
  clusterCut = 0.;
  maxMoveFlag = 1;
  if (transType_.compare("clustertrans") == 0) {
    maxMoveParam = 2e-3;
    targAcceptPer = 0.25;
    clusterCut = 4./3.;
  } else if (transType_.compare("clusterrotate") == 0) {
    maxMoveParam = 2e-3;
    targAcceptPer = 0.25;
    clusterCut = 4./3.;
  } else {
    ASSERT(0, "transformation type (" << transType_ << ") not recognized");
  }
}

void TrialCluster::attempt1_() {
  if (verbose_ == 1) {
    cout << std::setprecision(std::numeric_limits<double>::digits10+2)
         << "attempting " << transType_ << " " << pair_->peTot() << endl;
  }

  if (space()->nMol() <= 0) {
    trialMoveDecide_(0, 0);   // ensured rejection, however, criteria can update
    return void();
  }

  trialMoveRecordAll_(0);

  // compute and store cluster
  pair_->updateClusters(clusterCut);
  if (verbose_ == 1) cout << "number of clusters: " << space()->nClusters() << " n " << space()->nMol() << " cut " << clusterCut << endl;

  // rigidly translate or rotate clusters
  if (transType_.compare("clusterrotate") == 0) {
    space()->xClusterGen();
  }
  for (int i = 0; i < space()->nClusters(); ++i) {
    if (transType_.compare("clustertrans") == 0) {
      space()->randDispNoWrap(space()->clusterList()[i], maxMoveParam);
    } else if (transType_.compare("clusterrotate") == 0) {
      space()->randRotateMulti(space()->clusterList()[i], maxMoveParam,
                              pair_->sig());
    }
  }

  // wrap all molecules after transformation
  space()->wrapMol();

  // recompute and check cluster is same
  preFac_ = 1;
  vector<int> clusterMolOld = space()->clusterMol();
  if (space()->cellType() > 0) space()->updateCellofallMol();

  // update clusters, but do not accumulate stats if conf may be rejected
  space()->accumulateClusterVars_ = false;
  pair_->updateClusters(clusterCut);
  space()->accumulateClusterVars_ = true;

  for (int iMol = 0; iMol < space()->nMol(); ++iMol) {
    if (space()->clusterMol()[iMol] != clusterMolOld[iMol]) {
      preFac_ = 0.;
      if (verbose_ == 1) cout << "cluster changed" << endl;
    }
  }

  // compute forces and energy of new configuration
  if (preFac_ != 0) {
    de_ = pair_->allPartEnerForce(1) - peOld_;
    lnpMet_ = -criteria_->beta()*de_;
    reject_ = 0;
  } else {
    reject_ = 1;
    de_ = 0;
    lnpMet_ = std::numeric_limits<double>::min();
  }

  // accept or reject with bias prefactor
  if (criteria_->accept(lnpMet_, pair_->peTot() + de_, trialType_.c_str(),
                        reject_) == 1) {
    pair_->update(space()->listAtoms(), 0, "update");
    if (pair_->neighOn()) pair_->buildNeighList();
    trialAccept_();
  } else {
    space()->restoreAll();
    if (space()->cellType() > 0) space()->updateCellofallMol();
    // cout << "rejected " << transType_ << " " << de_ << endl;
    trialReject_();
  }
}

void TrialCluster::tuneParameters() {
  // determine limits and percentage changes
  double percent = 0.05;   // percentage change
  double upperLimit = 0;         // do not change above this limit
  double lowerLimit = 1e-5;
  if (transType_.compare("clustertrans") == 0) {
    upperLimit = space()->minl()/4.;
  } else if (transType_.compare("clusterrotate") == 0) {
    upperLimit = 1e1;
  } else {
    ASSERT(0, "unrecognized transType(" << transType_ << ")");
  }
  updateMaxMoveParam_(percent, upperLimit, lowerLimit, targAcceptPer);
}

shared_ptr<TrialCluster> makeTrialCluster(Pair *pair,
  Criteria *criteria, const char* transType) {
  return make_shared<TrialCluster>(pair, criteria, transType);
}

shared_ptr<TrialCluster> makeTrialCluster(const char* transType) {
  return make_shared<TrialCluster>(transType);
}

void clusterTrial(MC *mc, const char* type, const double clusterCut,
  const double maxMoveParam) {
  shared_ptr<TrialCluster> trial = make_shared<TrialCluster>(type);
  if (clusterCut != -1) trial->clusterCut = clusterCut;
  if (maxMoveParam != -1) trial->maxMoveParam = maxMoveParam;
  mc->initTrial(trial);
}
void clusterTrial(shared_ptr<MC> mc, const char* type, const double clusterCut,
  const double maxMoveParam) {
  clusterTrial(mc.get(), type, clusterCut, maxMoveParam);
}

}  // namespace feasst

