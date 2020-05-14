/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./trial_add.h"
#include "./mc.h"

namespace feasst {

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
  WARN(verbose_ == 1, "attempting to " << trialType_ << " " << molType_);
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

  double insertion_volume = 0.;
  if (confineFlag_ == 0) {
    // add molecule to entire box
    space()->addMol(molType_.c_str());
    insertion_volume = space()->volume();

    // obtain vector of particle IDs of inserted molecule
    mpart_ = space()->lastMolIDVec();

    if (twoParticle != 0) {
      space()->addMol(secondType.c_str());
      const vector<int> mpart2 = space()->lastMolIDVec();
      mpart_.insert(mpart_.end(), mpart2.begin(), mpart2.end());
    }
    //SWH: insert second molecule type twice
    if (threeParticle != 0) {
      space()->addMol(secondType.c_str());
      vector<int> mpart2 = space()->lastMolIDVec();
      mpart_.insert(mpart_.end(), mpart2.begin(), mpart2.end());
      space()->addMol(secondType.c_str());
      mpart2 = space()->lastMolIDVec();
      mpart_.insert(mpart_.end(), mpart2.begin(), mpart2.end());
    }
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
    insertion_volume = confine_volume_();
  }
  pair_->addPart();

  // if particle types are constrained to be equiMolar
  if (equiMolarCheckFails_(space()->mol().back(), 0)) {
//    cout << "constraint?" << endl;
    reject_ = 1;
  }

  if (avbOn_ && (reject_ != 1)) {
    if (space()->nMol() > 1) {
      ASSERT(confineFlag_ == 0, "AVB + confine may have issu with prefactor?");
      // HWH NOTE: Haven't implemented semigrand with AVB
      // select a random molecule that is different from mpart,
      // record its aggregation volume properties
      region_.assign("bonded");
      tmpart_ = space()->randMolDiff(space()->mol()[mpart_[0]]);
      vector<int> jneigh = pair_->neigh()[space()->mol()[tmpart_[0]]];
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
      preFac_ = vIn_*(nMol-1)/nMol/(nIn + 1);
      //preFac_ = vIn_*(space()->nMol()-1)/space()->nMol()/(nIn + 1);
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
    preFac_ = insertion_volume/static_cast<double>(nMolOfType);
//    cout << "prefac " << preFac_ << " vol " << insertion_volume << " n " << nMolOfType << endl;
    if (twoParticle != 0) {
      const int jMolIndex = space()->findAddMolListIndex(secondType);
      const int nMolOfjType = space()->nMolType()[jMolIndex];
      preFac_ *= insertion_volume/static_cast<double>(nMolOfjType);
//      cout << "prefac " << preFac_ << " vol " << insertion_volume << " n " << nMolOfjType << endl;
    }
    //SWH
    if (threeParticle != 0) {
      const int jMolIndex = space()->findAddMolListIndex(secondType);
      const int nMolOfjType = space()->nMolType()[jMolIndex];
      preFac_ *= insertion_volume*insertion_volume/(static_cast<double>(nMolOfjType)*(static_cast<double>(nMolOfjType)+1.0));
    }
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
//    cout << "lnpm " << lnpMet_ << endl;
    lnpMet_ += -criteria_->beta()*(de_ - def_)
            + log(criteria_->activ(iMolIndex));
//      cout << "lnpm " << lnpMet_ << " beta " << criteria_->beta() << " de " << de_ << " def " << def_ << " lnz " << log(criteria_->activ(iMolIndex)) << endl;
    if (twoParticle != 0) {
      const int jMolIndex = space()->findAddMolListIndex(secondType);
      lnpMet_ += log(criteria_->activ(jMolIndex));
//      cout << "lnpm " << lnpMet_ << " lnz " << log(criteria_->activ(jMolIndex)) << endl;
    }
    //SWH
    if (threeParticle != 0) {
      const int jMolIndex = space()->findAddMolListIndex(secondType);
      lnpMet_ += 2*log(criteria_->activ(jMolIndex));
    }
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
  int jMolIndex = 0;
  //SWH
  if (twoParticle != 0 || threeParticle != 0) {
    jMolIndex = space()->findAddMolListIndex(secondType);
  }
  if (header) {
    stat << "N" << iMolIndex << " ";
    if (twoParticle != 0 || threeParticle != 0) {
      stat << "N" << jMolIndex << " ";
    }
  } else {
    const int nMolOfType = space()->nMolType()[iMolIndex];
    stat << nMolOfType << " ";
    if (twoParticle != 0 || threeParticle != 0) {
      const int nMolOfjType = space()->nMolType()[jMolIndex];
      stat << nMolOfjType << " ";
    }
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

}  // namespace feasst
