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

namespace feasst {

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
  if (verbose_ == 1) std::cout << "attempting to " << trialType_ << " "
    << molType_ << std::endl;
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
      // select a random molecule tmpart as target for deletion, record its
      // aggregation volume properties
      region_.assign("bonded");
      tmpart_ = space()->randMol();
      vector<int> jneigh = pair_->neigh()[space()->mol()[tmpart_[0]]];
      //ASSERT(molType_.empty(), "haven't implemented semigrand for AVB");
      const int nMol = space()->nMol();
      if (!molType_.empty()) {
        //nMol = space()->nMolOfType(molType_);
        // remove neighbors which are not of molType_
        for (int i = static_cast<int>(jneigh.size()) - 1;
             i >= 0;
             --i) {
          if (space()->moltype()[jneigh[i]] != molType_) {
            jneigh.erase(jneigh.begin() + i);
          }
        }
      }
      const int nIn = static_cast<int>(jneigh.size());
      // cout << "nIn " << nIn << endl;
      if (nIn > 0) {
        preFac_ = nIn / vIn_ * nMol / (nMol - 1);
        //preFac_ = nIn / vIn_ * space()->nMol() / (space()->nMol() - 1);
        lnpMet_ = log(preFac_);
        mpart_ = space()->randMolSubset(jneigh);
      } else {
        reject_ = 1;
      }
    } else {
      // compute the deletion volume
      double deletion_volume = space()->volume();
      if (confineFlag_ != 0) {
        deletion_volume = confine_volume_();
      }

      // if no molType given, delete any molecule
      if (molType_.empty()) {
        preFac_ = static_cast<double>(space()->nMol())/deletion_volume;
        lnpMet_ = log(preFac_);
        mpart_ = space()->randMol();
        ASSERT(twoParticle == 0, "twoParticle not implemented without molType");
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
            int nMolOfjType = 0;
            if (twoParticle != 0) {
              const int jMolIndex = space()->findAddMolListIndex(secondType);
              nMolOfjType = space()->nMolType()[jMolIndex];
              if (nMolOfjType > 0) {
                const int jMol = space()->randMolofType(jMolIndex);
                if (jMol == -1) {
                  reject_ = 1;
                } else {
                  const vector<int> mpart2 = space()->imol2mpart(jMol);
                  mpart_.insert(mpart_.end(), mpart2.begin(), mpart2.end());
                  // sort mpart_ so that there are no issues if deleted
                  std::sort(mpart_.begin(), mpart_.end());
                }
              }
            }
            //SWH
            if (threeParticle != 0) {
              const int jMolIndex = space()->findAddMolListIndex(secondType);
              nMolOfjType = space()->nMolType()[jMolIndex];
              if (nMolOfjType > 0) {
                const int jMol = space()->randMolofType(jMolIndex);
                if (jMol == -1) {
                  reject_ = 1;
                } else {
                  vector<int> mpart2 = space()->imol2mpart(jMol);
                  mpart_.insert(mpart_.end(), mpart2.begin(), mpart2.end());
                  const int jMol2 = space()->randMolofType(jMolIndex);
                  if (jMol2 == -1 || jMol2 == jMol) {
                    reject_ = 1;
                  } else {
                    vector<int> mpart2 = space()->imol2mpart(jMol2);
                    mpart_.insert(mpart_.end(), mpart2.begin(), mpart2.end());
                    // sort mpart_ so that there are no issues if deleted
                    std::sort(mpart_.begin(), mpart_.end());
                  }
                }
              }
            }
            preFac_ = static_cast<double>(nMolOfType)/deletion_volume;
            if (twoParticle != 0) {
              preFac_ *= static_cast<double>(nMolOfjType)/deletion_volume;
            }
            //SWH
            if (threeParticle != 0) {
              preFac_ *= static_cast<double>(nMolOfjType)*(static_cast<double>(nMolOfjType)+1.0)/(deletion_volume*deletion_volume);
            }
            lnpMet_ = log(preFac_);
          }
        }
      }
    }

    // check if particle types are constrained to be equimolar
    if ( (reject_ != 1) &&
         (preFac_ != 0) &&
         (mpart_.size() > 0) ) {
      if (equiMolarCheckFails_(space()->mol()[mpart_[0]], 1)) {
        reject_ = 1;
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
      if (twoParticle != 0) {
        const int jMolIndex = space()->findAddMolListIndex(secondType);
        lnpMet_ -= log(criteria_->activ(jMolIndex));
      }
      //SWH
      if (threeParticle != 0) {
        const int jMolIndex = space()->findAddMolListIndex(secondType);
        lnpMet_ -= 2*log(criteria_->activ(jMolIndex));
      }
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

}  // namespace feasst

