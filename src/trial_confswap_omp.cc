/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#include "./trial_confswap_omp.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

TrialConfSwapOMP::TrialConfSwapOMP() : Trial() {
  defaultConstruction_();
}

TrialConfSwapOMP::TrialConfSwapOMP(
  Pair *pair,
  Criteria *criteria)
  : Trial(pair, criteria) {
  defaultConstruction_();
}

TrialConfSwapOMP::TrialConfSwapOMP(const char* fileName,
  Pair *pair,
  Criteria *criteria)
  : Trial(pair, criteria, fileName) {
  defaultConstruction_();
  string strtmp = fstos("orderTolerance", fileName);
  if (!strtmp.empty()) {
    orderTolerance_ = stod(strtmp);
  }
}

void TrialConfSwapOMP::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  file << "# orderTolerance " << orderTolerance_ << endl;
}

void TrialConfSwapOMP::defaultConstruction_() {
  className_.assign("TrialConfSwapOMP");
  trialType_.assign("move");
  verbose_ = 0;
  orderTolerance_ = 1e-4;
}

void TrialConfSwapOMP::attempt1_() {
  // obtain currentOrder, the current order parameter for the simulation
  double currentOrder = -1;
  if (orderType_.compare("nmol") == 0) {
    currentOrder = space()->nMol();
  } else if (orderType_.compare("pairOrder") == 0) {
    currentOrder = pair_->order();
  } else if (orderType_.compare("beta") == 0) {
    currentOrder = criteria_->beta();
  } else {
    ASSERT(0, "unrecognized macrostate type (" << orderType_ << ")");
  }

  // obtain procIndex, neighboring processors which overlap with the same
  // order parameter
  vector<int> procIndex;
  for (unsigned int i = 0; i < order_.size(); ++i) {
    const int index = order2index(currentOrder);
    if (index != -1) procIndex.push_back(index);
  }

  // if no processor overlaps, skip trial and don't count the attempt for
  // trial statistics
  if (procIndex.size() == 0) {
    --attempted_;
  } else {
    // randomly choose one of the processor indices to attempt a swap or store
    const int index = procIndex[uniformRanNum
      (0, static_cast<int>(procIndex.size()) - 1)];

    // OMP critical to prevent race condition between processors
    #pragma omp critical
    {
      // 1/2 chance to store configurations
      if (uniformRanNum() < 0.5) {
        // store current configuration but reject move
        confIntra_[index] = space()->cloneShrPtr();
        confIntra_[index]->cellOff();
        pe_[index] = pair_->peTot();
        --attempted_;
      // 1/2 chance to swap configurations on this processor with stored
      // neighboring proc
      } else {
        // check if inter proc configuration is stored
        TrialConfSwapOMP* trial = trialSwapInter_[index];
        const int indexInter = trial->order2index(currentOrder);
        ASSERT(indexInter != -1,
          "inter proc swap failed because indexInter == -1");
        Space *space_ptr = trial->confIntra()[indexInter].get();
        if ( (trial->confIntra()[indexInter] == NULL) ||
             (space()->natom() != space_ptr->natom()) ) {
          --attempted_;
        } else {
          const double peNew = trial->pe()[indexInter];
          de_ = peNew - pair_->peTot();
          lnpMet_ = space()->nMol()*dlnz_[index] - peNew*dbeta_[index]
                    - criteria_->beta()*de_;
          reject_ = 0;
          if (criteria_->accept(lnpMet_, pair_->peTot() + de_,
                                trialType_.c_str(), reject_) == 1) {
            trialAccept_();
            space()->swapPositions(space_ptr);
            if (space()->cellType() > 0) space()->buildCellList();
            if (pair_->neighOn()) pair_->buildNeighList();
            pair_->initEnergy();
          } else {
            trialReject_();
          }
        }
      }
    }
  }
}

void TrialConfSwapOMP::addProcOverlap(const double order,
  TrialConfSwapOMP* trial, const double dbeta, const double dlnz) {
  order_.push_back(order);
  pe_.push_back(0);
  trialSwapInter_.push_back(trial);
  confIntra_.resize(pe_.size());
  dbeta_.push_back(dbeta);
  dlnz_.push_back(dlnz);
}

int TrialConfSwapOMP::order2index(const double order) {
  int index = -1;
  for (unsigned int i = 0; i < order_.size(); ++i) {
    if (fabs(order - order_[i]) < orderTolerance_) index = i;
  }
  return index;
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_


