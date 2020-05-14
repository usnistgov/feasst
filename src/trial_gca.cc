/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./trial_gca.h"
#include "./mc.h"

namespace feasst {

TrialGCA::TrialGCA() : Trial() {
  defaultConstruction_();
}

TrialGCA::TrialGCA(
  Pair *pair,
  Criteria *criteria)
  : Trial(pair, criteria) {
  defaultConstruction_();
}

TrialGCA::TrialGCA(const char* fileName,
  Pair *pair,
  Criteria *criteria)
  : Trial(pair, criteria, fileName) {
  defaultConstruction_();
  maxMoveParam = fstod("maxMoveParam", fileName);
  targAcceptPer = fstod("targAcceptPer", fileName);
}

void TrialGCA::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  file << "# maxMoveParam " << maxMoveParam << endl;
  file << "# targAcceptPer " << targAcceptPer << endl;
}

void TrialGCA::defaultConstruction_() {
  className_.assign("TrialGCA");
  trialType_.assign("move");
  verbose_ = 0;
  maxMoveFlag = 1;
  maxMoveParam = 0.1;
  targAcceptPer = 10;
}

void TrialGCA::attempt1_() {
  if (verbose_ == 1) {
    cout << std::setprecision(std::numeric_limits<double>::digits10+2)
         << "attempting gca " << pair_->peTot() << endl;
  }
  if (space()->nMol() > 0) {
    attemptGCA_();
  } else {
    trialMoveDecide_(0, 0);   // ensured rejection, however, criteria can update
  }
}

void TrialGCA::tuneParameters() {
  // determine limits and percentage changes
  double percent = 0.05;   // percentage change
  double upperLimit = 0;         // do not change above this limit
  double lowerLimit = 1e-5;
  upperLimit = space()->minl()/2.;
  percent *= -1;
  updateMaxMoveParam_(percent, upperLimit, lowerLimit, targAcceptPer);
}

void TrialGCA::attemptGCA_() {
  if (space()->nMol() == 0) {
    reject_ = 1;
  } else {
    // // this line creates movie for each pivot in one GCA move
    // pair_->printXYZ("asdf",1);

    // turn on neighCutOn in pair class to record all neighbors
    // from multiPartEner
    pair_->initNeighCut(1);

    // turn on peMap for multiPartEner to record all pair-wise interactions,
    // in same order as neighCut
    pair_->initPEMap(1);

    // pivot a random molecule, and recursively select interacting molecules
    // to pivot
    const vector<int> mpart = space()->randMol();
    vector<int> pivotList;    // list of molecules that are pivoted
    vector<int> rejectNeigh;  // list of molecules that rejected a pivot
    vector<double> rejectDe;  // list of energy of interaction of rejected pivot

    // position of the pivot
    const vector<double> xPivot = space()->randPosition(
      space()->mol()[mpart[0]], maxMoveParam);
    recursivePivot_(mpart, xPivot, &pivotList, &rejectNeigh, &rejectDe);

    // turn off neighCutOn in pair class to record all neighbors from
    // multiPartEner
    pair_->initNeighCut(0);

    // turn off peMap
    pair_->initPEMap(0);

    // compute change in energy
    // cout << "number of rejects " << int(rejectNeigh.size()) << endl;
    // cout << "number pivoted " << int(pivotList.size()) << endl;
    for (unsigned int iNeigh = 0; iNeigh < rejectNeigh.size(); ++iNeigh) {
      if (!findInList(rejectNeigh[iNeigh], pivotList)) {
        // cout << "derej " << rejectNeigh[iNeigh] << " "
        //      << rejectDe[iNeigh] << endl;
        de_ += rejectDe[iNeigh];
      }
    }
    // cout << "change of energ " << de_ << endl;

    // always accept the trial
    pair_->update(de_);
  }
  if (criteria_->accept(0, pair_->peTot(), "move", reject_) == 1) {
    trialAccept_();
  } else {
    trialReject_();
  }
}

/**
 * recursively pivot molecule iMol
 */
void TrialGCA::recursivePivot_(
  const vector<int> mpart,        //!< list of sites of pivoting molecule
  const vector<double> &xPivot,   //!< position of pivot
  vector<int> *pivotList,    //!< list of molecules that have been pivoted
  vector<int> *rejectNeigh,  //!< molecules that were rejected from pivot
  vector<double> *rejectDe   //!< interactions from rejected pivots
  ) {
  // check if molecule is already in the pivot list
  const int iMol = space()->mol()[mpart.front()];
  if (!findInList(iMol, *pivotList)) {
    // cout << "attempting to pivot iMol " << iMol << endl;

    // count number of pivoted molecules for MC class tuning of maxMoveParam
    if (pivotList->size() > 0) ++accepted_;

    // add molecule to the pivot list
    pivotList->push_back(iMol);

    // find all interacting neighbors of iMol at location before pivot,
    // and their interactions
    peOld_ = pair_->multiPartEner(mpart, 0);
    pair_->update(mpart, 4, "store");

    // list of interacting neighbors in old and new config
    vector<int> neighOld, neigh;

    // list of energy interactions in same order as neigh
    vector<double> peMapOld, peMap;
    pair_->neighCutMolPEMap(neighOld, peMapOld);

    // reflect iMol about a point
    space()->pointReflectMol(iMol, xPivot);
    space()->wrap(mpart);
    if (space()->cellType() > 0) space()->updateCellofiMol(iMol);

    // // this line creates movie for each pivot in one GCA move
    //  pair_->printXYZ("asdf",0);

    // find all interacting neighbors of iMol at new location,
    // and their interactions
    pair_->multiPartEner(mpart, 0);
    pair_->update(mpart, 4, "update");
    pair_->neighCutMolPEMap(neigh, peMap);

    // generate list of molecules that will be considered for a recursive
    // pivot, and their next ixn difference
    vector<int> newPivots = neighOld;
    vector<double> newPEMap = peMapOld;
    for (unsigned int i = 0; i < newPEMap.size(); ++i) {
      newPEMap[i] *= -1;
    }

    // combine set of neigh and neighOld without double counting,
    // and compute interaction changes
    //  do this by looping through neigh, and not adding those which have
    //  already been counted
    for (unsigned int ii = 0; ii < neigh.size(); ++ii) {
      const int iNeigh = neigh[ii];
      int index = -1;
      if (findInList(iNeigh, neighOld, &index)) {
        newPEMap[index] += peMap[ii];
      } else {
        newPivots.push_back(iNeigh);
        newPEMap.push_back(peMap[ii]);
      }
    }

    // use interaction energy for probability of adding to pivot
    for (unsigned int ii = 0; ii < newPivots.size(); ++ii) {
      const int jMol = newPivots[ii];
      if (uniformRanNum() < 1 - exp(-criteria_->beta()*newPEMap[ii])) {
        recursivePivot_(space()->imol2mpart(jMol), xPivot, pivotList,
                       rejectNeigh, rejectDe);
      } else {
        rejectNeigh->push_back(jMol);
        rejectDe->push_back(newPEMap[ii]);
      }
    }
  }
}

shared_ptr<TrialGCA> makeTrialGCA(Pair *pair,
  Criteria *criteria) {
  return make_shared<TrialGCA>(pair, criteria);
}

shared_ptr<TrialGCA> makeTrialGCA() {
  return make_shared<TrialGCA>();
}

void gcaTrial(MC *mc, const int nMolTarg) {
  shared_ptr<TrialGCA> trial = make_shared<TrialGCA>();
  if (nMolTarg != -1) trial->targAcceptPer = nMolTarg;
  mc->initTrial(trial);
}

void gcaTrial(shared_ptr<MC> mc, const int nMolTarg) {
  gcaTrial(mc.get(), nMolTarg);
}

}  // namespace feasst



