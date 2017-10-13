/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#include "./trial_add.h"
#include "./mc.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

TrialAdd::TrialAdd(const char* molType)
  : Trial(),
    molType_(molType) {
  defaultConstruction_();
}

TrialAdd::TrialAdd(
  Pair *pair,
  Criteria *criteria,
  const char* molType)
  : Trial(pair, criteria),
    molType_(molType) {
  defaultConstruction_();
}

TrialAdd::TrialAdd(const char* fileName,
  Pair *pair,
  Criteria *criteria)
  : Trial(pair, criteria, fileName) {
  defaultConstruction_();
  molType_ = fstos("molType", fileName);
}

void TrialAdd::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  file << "# molType " << molType_ << endl;
}

void TrialAdd::defaultConstruction_() {
  className_.assign("TrialAdd");
  trialType_.assign("add");
  verbose_ = 0;
  molid_ = -1;
}

TrialAdd* TrialAdd::clone(Pair *pair, Criteria *criteria) const {
  TrialAdd* t = new TrialAdd(*this);
  t->reconstruct(pair, criteria);
  return t;
}
shared_ptr<TrialAdd> TrialAdd::cloneShrPtr(
  Pair* pair, Criteria* criteria) const {
  return(std::static_pointer_cast<TrialAdd, Trial>
    (cloneImpl(pair, criteria)));
}
shared_ptr<Trial> TrialAdd::cloneImpl(
  Pair *pair, Criteria *criteria) const {
  shared_ptr<TrialAdd> t = make_shared<TrialAdd>(*this);
  t->reconstruct(pair, criteria);
  return t;
}


void TrialAdd::attempt1_() {
  WARN(verbose_ == 1, "attempting to " << trialType_);
  ASSERT((pair_->atomCut() != 1) || (space()->nMol() == space()->natom()) ||
         (!avbOn_), "this class assumes atomCut(" << pair_->atomCut()
         << ") == 0 when avb is on");

  // initialize molid_ if not already initialized from default value
  if (molid_ == -1) {
    molid_ = space()->findAddMolListIndex(molType_);
    ASSERT(molid_ != -1, "molType(" << molType_ << " not initialized.");
    std::stringstream ss;
    ss << "add" << molid_;
    trialType_.assign(ss.str());
  }

  if (confineFlag_ == 0) {
    // add molecule to entire box
    space()->addMol(molType_.c_str());

    // obtain vector of particle IDs of inserted molecule
    mpart_ = space()->lastMolIDVec();
  } else {
    // otherwise, add molecule only to confines
    bool inRegion = false;
    int tries = 0, maxTries = 1e4;
    while (!inRegion && (tries < maxTries)) {
      space()->addMol(molType_.c_str());
      mpart_ = space()->lastMolIDVec();
      if ( (space()->x(mpart_[0], confineDim_) <= confineUpper_) &&
           (space()->x(mpart_[0], confineDim_) >= confineLower_) ) {
        inRegion = true;
      } else {
        space()->delPart(mpart_);
      }
      ++tries;
    }
    ASSERT(tries != maxTries, "maximum number of attempts in confine");
  }
  pair_->addPart();

  // if particle types are constrained to be equiMolar
  if (space()->equiMolar() >= 1) {
    const int iMol = space()->mol().back(),
      iMolType = space()->molid()[iMol];
    if (space()->equiMolar() == 1) {
      const int nimt = space()->nMolType()[iMolType];

      // check if a different moltype has 2+ less in number
      for (int imt = 0; imt < space()->nMolTypes(); ++imt) {
        if (imt != iMolType) {
          if (nimt > space()->nMolType()[imt]+1) reject_ = 1;
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

  if (avbOn_ && (reject_ != 1)) {
    if (space()->nMol() > 1) {
      // HWH NOTE: Haven't implemented semigrand with AVB
      // select a random molecule that is different from mpart,
      // record its aggregation volume properties
      region_.assign("bonded");
      tmpart_ = space()->randMolDiff(space()->mol()[mpart_[0]]);
      const int nIn = static_cast<int>(pair_->neigh()
        [space()->mol()[tmpart_[0]]].size());
      preFac_ = vIn_*(space()->nMol()-1)/space()->nMol()/(nIn + 1);
      lnpMet_ = log(preFac_);

      // move mpart to bonded region of tmpart
      space()->avb(mpart_, tmpart_, rAbove_, rBelow_, region_.c_str());
      if (space()->cellType() > 0) {
        space()->updateCellofiMol(space()->mol()[mpart_.front()]);
      }
    }
  } else {
    const int iMolIndex = space()->findAddMolListIndex(molType_);
    const int nMolOfType = space()->nMolType()[iMolIndex];
    preFac_ = space()->vol()/static_cast<double>(nMolOfType);
    lnpMet_ = log(preFac_);
  }
  space()->wrap(mpart_);

  // multiple first bead insertion modifies def_, preFac_, tmpart_ for avb
  if ( (nf_ > 1) && (preFac_ != 0) && (reject_ != 1) ) {
    const double w = multiFirstBead_(3);
    lnpMet_ += log(w);
    preFac_ *= w;
    if (space()->cellType() > 0) {
      space()->updateCellofiMol(space()->mol()[mpart_.front()]);
    }
  }

  // record energy contribution of selected particle
  if ( (preFac_ != 0) && (reject_ != 1) ) {
    de_ = pair_->multiPartEner(mpart_, 3);
    pair_->update(mpart_, 3, "store");
    const int iMolIndex = space()->findAddMolListIndex(molType_);
    lnpMet_ += -criteria_->beta()*(de_ - def_)
            + log(criteria_->activ(iMolIndex));
    reject_ = 0;
  } else {
    reject_ = 1;
    de_ = 0;
    lnpMet_ = std::numeric_limits<double>::min();
  }

  if (criteria_->accept(lnpMet_, pair_->peTot() + de_, trialType_.c_str(),
                        reject_) == 1) {
    trialAccept_();
    pair_->update(mpart_, 3, "update");
    WARN(verbose_ == 1, "insertion accepted " << de_);

  // if not accepted, remove molecule, assuming a molecule is described
  // by sequentially listed particles
  } else {
    pair_->delPart(mpart_);
    space()->delPart(mpart_);
    WARN(verbose_ == 1, "insertion rejected " << de_);
    trialReject_();
  }
}

string TrialAdd::printStat(const bool header) {
  stringstream stat;
  stat << Trial::printStat(header);
  const int iMolIndex = space()->findAddMolListIndex(molType_);
  if (header) {
    stat << "N" << iMolIndex << " ";
  } else {
    const int nMolOfType = space()->nMolType()[iMolIndex];
    stat << nMolOfType << " ";
  }
  return stat.str();
}

shared_ptr<TrialAdd> makeTrialAdd(Pair *pair, Criteria *criteria,
  const char* molType) {
  return make_shared<TrialAdd>(pair, criteria, molType);
}

shared_ptr<TrialAdd> makeTrialAdd(const char* molType) {
  return make_shared<TrialAdd>(molType);
}

void addTrial(MC *mc, const char* moltype) {
  shared_ptr<TrialAdd> trial = make_shared<TrialAdd>(moltype);
  mc->initTrial(trial);
}

void addTrial(shared_ptr<MC> mc, const char* moltype) {
  addTrial(mc.get(), moltype);
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_
