#include "./trial_delete.h"
#include <stdio.h>
#include <iostream>
#include <numeric>
#include "./functions.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

TrialDelete::TrialDelete() : Trial() {
  defaultConstruction();
  molType_.assign("");
}
TrialDelete::TrialDelete(const char* molType)
  : Trial(),
  molType_(molType) {
  defaultConstruction();
}
TrialDelete::TrialDelete(Space *space,
  Pair *pair,
  Criteria *criteria)
  : Trial(space, pair, criteria) {
  defaultConstruction();
  molType_.assign("");
}
TrialDelete::TrialDelete(Space *space,
  Pair *pair,
  Criteria *criteria,
  const char* molType)
  : Trial(space, pair, criteria),
  molType_(molType) {
  defaultConstruction();
}
TrialDelete::TrialDelete(const char* fileName,
  Space *space,
  Pair *pair,
  Criteria *criteria)
  : Trial(space, pair, criteria, fileName) {
  defaultConstruction();
  molType_ = fstos("molType", fileName);
}

/**
 * default constructor
 */
void TrialDelete::defaultConstruction() {
  className_.assign("TrialDelete");
  trialType_.assign("del");
  molid_ = -1;
  verbose_ = 0;
}

/**
 * write restart file
 */
void TrialDelete::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  if (!molType_.empty()) {
    file << "# molType " << molType_ << endl;
  }
}

/**
 * Attempt deletion of a molecule.
 */
void TrialDelete::attempt1() {
  ASSERT((pair_->atomCut() != 1) || (space_->nMol() == space_->natom()) ||
         (!avbOn_), "this class assumes atomCut(" << pair_->atomCut()
         << ") == 0 when avb is on");
  if (verbose_ == 1) std::cout << "attempting to " << trialType_ << std::endl;

  // initialize molid_ if not already initialized from default value
  if (molid_ == -1) {
    if (molType_.empty()) {
      molid_ = 0;
    } else {
      molid_ = space_->findAddMolListIndex(molType_);
    }
    ASSERT(molid_ != -1, "molType(" << molType_ << " not initialized.");
    std::stringstream ss;
    ss << "del" << molid_;
    trialType_.assign(ss.str());
  }

  if (space_->nMol() <= 0) {
    if (verbose_ == 1) {
      cout << "deletion rejected because no molecules " << de_ << endl;
    }
  } else if ( (avbOn_) && (space_->nMol() <= 1) ) {
    if (verbose_ == 1) {
      cout << "deletion rejected because not enough molecules for avb " << de_
           << std::endl;
    }
  } else {
    // select a molecule to delete, mpart
    if (avbOn_) {
      ASSERT(molType_.empty(), "haven't implemented semigrand for AVB");
      // select a random molecule tmpart as target for deletion, record its
      // aggregation volume properties
      region_.assign("bonded");
      tmpart_ = space_->randMol();
      vector<int> jneigh = pair_->neigh()[space_->mol()[tmpart_[0]]];
      const int nIn = static_cast<int>(jneigh.size());
      if (nIn > 0) {
        preFac_ = nIn / vIn_ * space_->nMol() / (space_->nMol() - 1);
        lnpMet_ = log(preFac_);
        mpart_ = space_->randMolSubset(jneigh);
      }
    } else {
      // if no molType given, delete any molecule
      if (molType_.empty()) {
        preFac_ = static_cast<double>(space_->nMol())/space_->vol();
        lnpMet_ = log(preFac_);
        mpart_ = space_->randMol();
    
      // otherwise, delete only molType
      } else {
        const int iMolIndex = space_->findAddMolListIndex(molType_);
        const int nMolOfType = space_->nMolType()[iMolIndex];
        if (nMolOfType > 0) {
          const int iMol = space_->randMolofType(iMolIndex);
          if (iMol == -1) {
            reject_ = 1;
          } else {
            mpart_ = space_->imol2mpart(iMol);
            preFac_ = static_cast<double>(nMolOfType)/space_->vol();
            lnpMet_ = log(preFac_);
          }
        }
      }
    }

    // check if particle types are constrained to be equimolar
    if (space_->equiMolar() >= 1) {
      const int iMol = space_->mol()[mpart_[0]],
        iMolType = space_->molid()[iMol];
      if (space_->equiMolar() == 1) {
        const int nimt = space_->nMolType()[iMolType];

        // check if any case where there is already a moltype with more
        // molecules
        for (int imt = 0; imt < space_->nMolTypes(); ++imt) {
          if (imt != iMolType) {
            if (nimt < space_->nMolType()[imt]) reject_ = 1;
          }
        }
      } else if (space_->equiMolar() == 2) {
        if (space_->nMol() % 2 == 0) {
          if (iMolType == 0) reject_ = 1;
        }
        if (space_->nMol() % 2 == 1) {
          if (iMolType == 1) reject_ = 1;
        }
      } else if (space_->equiMolar() == 3) {
        if (space_->nMol() % 2 == 0) {
          if (iMolType == 1) reject_ = 1;
        }
        if (space_->nMol() % 2 == 1) {
          if (iMolType == 0) reject_ = 1;
        }
      }
    }

    // catch for nIn == 0, mpart_ not defined
    if (preFac_ != 0) {
      // multiple first beads modifies def_, preFac_, tmpart for avb
      if (nf_ > 1) {
        const double w = multiFirstBead(2);
        lnpMet_ += log(w);
        preFac_ *= w;
      }

      // record energy contribution of molecule to delete
      de_ = -1. * pair_->multiPartEner(mpart_, 2);
      pair_->update(mpart_, 2, "store");
      int iMolIndex = -1;
      if (molType_.empty()) {
        iMolIndex = 0;
      } else {
        iMolIndex = space_->findAddMolListIndex(molType_);
      }
      lnpMet_ += -criteria_->beta()*(de_ + def_)
              - log(criteria_->activ(iMolIndex));
    }
  }

  if (preFac_ == 0) {
    reject_ = 1;
    de_ = 0;
    lnpMet_ = std::numeric_limits<double>::min();
  }

  // acceptance criteria
  if (criteria_->accept(lnpMet_, pair_->peTot() + de_,
                        trialType_.c_str(), reject_) == 1) {
    // remove molecule, assuming a molecule is described by sequentially listed
    // particles
    pair_->delPart(mpart_);
    space_->delPart(mpart_);
    pair_->update(mpart_, 2, "update");
    trialAccept();
    if (verbose_ == 1) cout << "deletion accepted " << de_ << std::endl;

  // if not accepted, restore
  } else {
    if (verbose_ == 1) std::cout << "deletion rejected " << de_ << std::endl;
    trialReject();
  }
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_



