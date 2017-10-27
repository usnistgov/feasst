/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./trial_delete.h"
#include <stdio.h>
#include <iostream>
#include <numeric>
#include "./functions.h"
#include "./mc.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

TrialDelete::TrialDelete() : Trial() {
  defaultConstruction_();
  molType_.assign("");
}

TrialDelete::TrialDelete(const char* molType)
  : Trial(),
  molType_(molType) {
  defaultConstruction_();
}

TrialDelete::TrialDelete(
  Pair *pair,
  Criteria *criteria)
  : Trial(pair, criteria) {
  defaultConstruction_();
  molType_.assign("");
}

TrialDelete::TrialDelete(
  Pair *pair,
  Criteria *criteria,
  const char* molType)
  : Trial(pair, criteria),
  molType_(molType) {
  defaultConstruction_();
}

TrialDelete::TrialDelete(const char* fileName,
  Pair *pair,
  Criteria *criteria)
  : Trial(pair, criteria, fileName) {
  defaultConstruction_();
  molType_ = fstos("molType", fileName);
}

void TrialDelete::defaultConstruction_() {
  className_.assign("TrialDelete");
  trialType_.assign("del");
  molid_ = -1;
  verbose_ = 0;
}

void TrialDelete::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  if (!molType_.empty()) {
    file << "# molType " << molType_ << endl;
  }
}

void TrialDelete::attempt1_() {
  ASSERT((pair_->atomCut() != 1) || (space()->nMol() == space()->natom()) ||
         (!avbOn_), "this class assumes atomCut(" << pair_->atomCut()
         << ") == 0 when avb is on");
  if (verbose_ == 1) std::cout << "attempting to " << trialType_ << std::endl;

  // initialize molid_ if not already initialized from default value
  if (molid_ == -1) {
    if (molType_.empty()) {
      molid_ = 0;
    } else {
      molid_ = space()->findAddMolListIndex(molType_);
    }
    ASSERT(molid_ != -1, "molType(" << molType_ << " not initialized.");
    std::stringstream ss;
    ss << "del" << molid_;
    trialType_.assign(ss.str());
  }

  if (space()->nMol() <= 0) {
    if (verbose_ == 1) {
      cout << "deletion rejected because no molecules " << de_ << endl;
    }
  } else if ( (avbOn_) && (space()->nMol() <= 1) ) {
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
      tmpart_ = space()->randMol();
      vector<int> jneigh = pair_->neigh()[space()->mol()[tmpart_[0]]];
      const int nIn = static_cast<int>(jneigh.size());
      if (nIn > 0) {
        preFac_ = nIn / vIn_ * space()->nMol() / (space()->nMol() - 1);
        lnpMet_ = log(preFac_);
        mpart_ = space()->randMolSubset(jneigh);
      }
    } else {
      // if no molType given, delete any molecule
      if (molType_.empty()) {
        preFac_ = static_cast<double>(space()->nMol())/space()->vol();
        lnpMet_ = log(preFac_);
        mpart_ = space()->randMol();

      // otherwise, delete only molType
      } else {
        const int iMolIndex = space()->findAddMolListIndex(molType_);
        const int nMolOfType = space()->nMolType()[iMolIndex];
        if (nMolOfType > 0) {
          const int iMol = space()->randMolofType(iMolIndex);
          if (iMol == -1) {
            reject_ = 1;
          } else {
            mpart_ = space()->imol2mpart(iMol);
            preFac_ = static_cast<double>(nMolOfType)/space()->vol();
            lnpMet_ = log(preFac_);
          }
        }
      }
    }

    // check if particle types are constrained to be equimolar
    if (space()->equiMolar() >= 1) {
      const int iMol = space()->mol()[mpart_[0]],
        iMolType = space()->molid()[iMol];
      if (space()->equiMolar() == 1) {
        const int nimt = space()->nMolType()[iMolType];

        // check if any case where there is already a moltype with more
        // molecules
        for (int imt = 0; imt < space()->nMolTypes(); ++imt) {
          if (imt != iMolType) {
            if (nimt < space()->nMolType()[imt]) reject_ = 1;
          }
        }
      } else if (space()->equiMolar() == 2) {
        if (space()->nMol() % 2 == 0) {
          if (iMolType == 0) reject_ = 1;
        }
        if (space()->nMol() % 2 == 1) {
          if (iMolType == 1) reject_ = 1;
        }
      } else if (space()->equiMolar() == 3) {
        if (space()->nMol() % 2 == 0) {
          if (iMolType == 1) reject_ = 1;
        }
        if (space()->nMol() % 2 == 1) {
          if (iMolType == 0) reject_ = 1;
        }
      }
    }

    // check if deletion is confined to a region
    if ( (confineFlag_ == 1) && (preFac_ != 0) ) {
      if ( (space()->x(mpart_[0], confineDim_) > confineUpper_) ||
           (space()->x(mpart_[0], confineDim_) < confineLower_) ) {
        preFac_ = 0.;
        reject_ = 1;
      }
    }

    // catch for nIn == 0, mpart_ not defined
    if (preFac_ != 0) {
      // multiple first beads modifies def_, preFac_, tmpart for avb
      if (nf_ > 1) {
        const double w = multiFirstBead_(2);
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
        iMolIndex = space()->findAddMolListIndex(molType_);
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
    space()->delPart(mpart_);
    pair_->update(mpart_, 2, "update");
    trialAccept_();
    if (verbose_ == 1) cout << "deletion accepted " << de_ << std::endl;

  // if not accepted, restore
  } else {
    if (verbose_ == 1) std::cout << "deletion rejected " << de_ << std::endl;
    trialReject_();
  }
}

shared_ptr<TrialDelete> makeTrialDelete(Pair *pair,
  Criteria *criteria, const char* molType) {
  return make_shared<TrialDelete>(pair, criteria, molType);
}

shared_ptr<TrialDelete> makeTrialDelete(const char* molType) {
  return make_shared<TrialDelete>(molType);
}

shared_ptr<TrialDelete> makeTrialDelete(Pair *pair,
  Criteria *criteria) {
  return make_shared<TrialDelete>(pair, criteria);
}

shared_ptr<TrialDelete> makeTrialDelete() {
  return make_shared<TrialDelete>();
}

void deleteTrial(MC *mc, const char* moltype) {
  shared_ptr<TrialDelete> trial = make_shared<TrialDelete>(moltype);
  mc->initTrial(trial);
}

void deleteTrial(shared_ptr<MC> mc, const char* moltype) {
  deleteTrial(mc.get(), moltype);
}

void deleteTrial(MC *mc) {
  shared_ptr<TrialDelete> trial = make_shared<TrialDelete>();
  mc->initTrial(trial);
}

void deleteTrial(shared_ptr<MC> mc) {
  deleteTrial(mc.get());
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

