#include "./trial_configBias.h"
#include "./mc.h"

namespace feasst {

TrialConfigBias::TrialConfigBias(const char* molType)
  : Trial(),
    molType_(molType) {
}
TrialConfigBias::TrialConfigBias(Space *space,
  Pair *pair,
  Criteria *criteria,
  const char* molType)   //!< type of molecule to add
  : Trial(space, pair, criteria),
    molType_(molType) {
  // // this line creates a movie
  // pair_->printxyz("asdf", 1);
}
TrialConfigBias::TrialConfigBias(const char* fileName,
  Space *space,
  Pair *pair,
  Criteria *criteria)
  : Trial(space, pair, criteria, fileName) {
}

/**
 * clone design pattern
 */
TrialConfigBias* TrialConfigBias::clone
  (Space* space, Pair *pair, Criteria *criteria) const {
  TrialConfigBias* t = new TrialConfigBias(*this);
  t->reconstruct(space, pair, criteria);
  return t;
}
shared_ptr<TrialConfigBias> TrialConfigBias::cloneShrPtr
  (Space* space, Pair* pair, Criteria* criteria) const {
  return(std::static_pointer_cast<TrialConfigBias, Trial>
    (cloneImpl(space, pair, criteria)));
}
shared_ptr<Trial> TrialConfigBias::cloneImpl
  (Space* space, Pair *pair, Criteria *criteria) const {
  shared_ptr<TrialConfigBias> t = make_shared<TrialConfigBias>(*this);
  t->reconstruct(space, pair, criteria);
  return t;
}

/**
 * add generic growth step
 */
void TrialConfigBias::addStepGeneric(const int iAtom, const int n) {
  stepAtom_.push_back(iAtom);
  stepTrials_.push_back(n);
  stepParam_.resize(stepAtom_.size());
  stepList_.resize(stepAtom_.size());
}

/**
 * add growth step, placing iAtom randomly in box n times
 */
void TrialConfigBias::addStepBox(const int iAtom, const int n) {
  addStepGeneric(iAtom, n);
  string name("box");
  stepName_.push_back(name);
}

/**
 * update stepParam, using current values from space class
 *  also, take into account step order
 */
void TrialConfigBias::updateStepParams() {
  stepParam_.resize(nSteps());
  shared_ptr<Space> s = space_->findAddMolInList(molType_);
  for (int iStep = 0; iStep < nSteps(); ++iStep) {
    const int iAtom = stepAtom_[iStep];
    int jAtom = -1, kAtom = -1, lAtom = -1;
    if (static_cast<int>(stepList_[iStep].size() > 0))  {
      jAtom = stepList_[iStep][0];
    }
    if (static_cast<int>(stepList_[iStep].size() > 1)) {
      kAtom = stepList_[iStep][1];
    }
    if (static_cast<int>(stepList_[iStep].size() > 2)) {
      lAtom = stepList_[iStep][2];
    }
    if (stepName_[iStep].compare("bond") == 0) {
      stepParam_[iStep] = s->bondParams(iAtom, jAtom);

    } else if (stepName_[iStep].compare("angle") == 0) {
      vector<double> bParams, aParams, allParams;
      bParams = s->bondParams(iAtom, jAtom);
      for (unsigned int i = 0; i < bParams.size(); ++i) {
        allParams.push_back(bParams[i]);
      }
      aParams = s->angleParams(iAtom, jAtom, kAtom);
      for (unsigned int i = 0; i < aParams.size(); ++i) {
        allParams.push_back(aParams[i]);
      }
      stepParam_[iStep] = allParams;

    } else if (stepName_[iStep].compare("branch") == 0) {
      vector<double> bParams, aParams, allParams;
      bParams = s->bondParams(jAtom, kAtom);
      for (unsigned int i = 0; i < bParams.size(); ++i) {
        allParams.push_back(bParams[i]);
      }
      bParams = s->bondParams(iAtom, kAtom);
      for (unsigned int i = 0; i < bParams.size(); ++i) {
        allParams.push_back(bParams[i]);
      }
      aParams = s->angleParams(iAtom, kAtom, lAtom);
      for (unsigned int i = 0; i < aParams.size(); ++i) {
        allParams.push_back(aParams[i]);
      }
      aParams = s->angleParams(jAtom, kAtom, lAtom);
      for (unsigned int i = 0; i < aParams.size(); ++i) {
        allParams.push_back(aParams[i]);
      }
      aParams = s->angleParams(iAtom, kAtom, jAtom);
      for (unsigned int i = 0; i < aParams.size(); ++i) {
        allParams.push_back(aParams[i]);
      }
      stepParam_[iStep] = allParams;

    } else if ( (stepName_[iStep].compare("box") != 0) &&
                (stepName_[iStep].compare("avb") != 0) &&
                (stepName_[iStep].compare("freeCOM") != 0) ) {
      ASSERT(0, "unrecognized stepName(" << stepName_[iStep]
             << ") in updateStepParams");
    }
  }

  // find steps with 1 trial and make them consecutive with the previous
  stepOrder_.resize(nSteps(), 1);
  for (int iStep = nSteps() - 1; iStep > 0; --iStep) {
    // if iStep has 1 trial move, remove step order here and add it to iStep-1
    if (stepTrials_[iStep] == 1) {
      stepOrder_[iStep-1] += stepOrder_[iStep];
      stepOrder_[iStep] = 0;
    }
  }
}

/**
 * add growth step, placing iAtom randomly along bond with jAtom
 */
void TrialConfigBias::addStepBond(
  const int iAtom, const int jAtom, const int n) {
  addStepGeneric(iAtom, n);
  string name("bond");
  stepName_.push_back(name);

  vector<int> list(1); list[0] = jAtom;
  stepList_.back() = list;
}

/**
 * add growth step, placing iAtom randomly along angle with jAtom and kAtom
 */
void TrialConfigBias::addStepAngle(
  const int iAtom, const int jAtom, const int kAtom, const int n) {
  addStepGeneric(iAtom, n);
  string name("angle");
  stepName_.push_back(name);

  vector<int> list(2); list[0] = jAtom; list[1] = kAtom;
  stepList_.back() = list;
}

/**
 * add growth step, placing iAtom and jAtom along branch with kAtom and lAtom, ntimes
 *  first, jAtom is placed using angle with kAtom and lAtom
 *  then, iAtom is placed according to both angles ikl and ikj
 *
 *               l
 *               .
 *       (tjkl)  .  (tikl)
 *               k
 *             .   .
 *          .        .
 *       (j)   (tikj) (i)
 *
 * when comparing with Space::setAtomInBranch
 *  a1 = lAtom
 *  a2 = jAtom
 *  a3 = iAtom
 *  a4 = kAtom
 */
void TrialConfigBias::addStepBranch(const int iAtom, const int jAtom,
  const int kAtom, const int lAtom, const int n) {
  addStepGeneric(iAtom, n);
  string name("branch");
  stepName_.push_back(name);

  vector<int> list(3); list[0] = jAtom; list[1] = kAtom; list[2] = lAtom;
  stepList_.back() = list;
}

/**
 * add growth step, placing iAtom at Center-of-Mass for free
 */
void TrialConfigBias::addStepFreeCOM(const int iAtom) {
  addStepGeneric(iAtom, 1);
  string name("freeCOM");
  stepName_.push_back(name);
}

/**
 * add aggregation volume bias (avb) step, placing iAtom in the avb of the list of target atoms, tmpart, with avb defined by rAbove and rBelow, n times
 */
void TrialConfigBias::addStepAVB(const int iAtom, const vector<int> tmpart,
  const double rAbove, const double rBelow, const int n) {
  addStepGeneric(iAtom, n);
  string name("avb");
  stepName_.push_back(name);
  stepList_.back() = tmpart;
  initAVB(rAbove, rBelow);
}

/**
 * Attempt to regrow/build molecule
 */
void TrialConfigBias::attempt1() {
  // upon every attempt, check to see if stepParam needs to be updated
  //  due to change in order parameter of pair class
  if (pair_->orderName().compare("angle0") == 0) {
    const double pairOrder = pair_->order();
    if (pairOrder != pairOrderOld_) {
      updateStepParams();
      pairOrderOld_ = pairOrder;
    }
  }

  if (growType.compare("regrow") == 0) {
    trialType_.assign("move");
    mout_("verbose", std::ostringstream().flush() << "attempting to "
          << trialType_);

    double lnwo, eo, lnwn, en, eof = 0;
    if (space_->nMol() == 0) {
      reject_ = 1;
    } else {
      // select a random molecule to regrow
      mpart_ = mpartall_ = space_->randMol();    // select a random molecule
      mpartCBparse();      // reduce mpart_ to only particles involved in CB

      // store old configuration
      // For AVB, buildAllStep may lead to a change in mpart_ and mpartall_
      buildAllStep(eo, lnwo, 0);
      if (dualCut_ > 0) {
        eof = pair_->multiPartEner(mpart_, 0);
        pair_->update(mpart_, 0, "store");
      }
      space_->xStore(mpart_);

      // compute new rosenbluth weights and new external energy
      buildAllStep(en, lnwn, 1);
    }

    // acceptance criteria
    if (reject_ == 1) {
      de_ = 0;
      lnpMet_ = std::numeric_limits<double>::min();
    } else {
      de_ = en - eo;
      lnpMet_ = lnwn - lnwo;

      // if dual cut config bias, compute full energy term
      if (dualCut_ > 0) {
        const double enf = pair_->multiPartEner(mpart_, 1);
        lnpMet_ -= criteria_->beta()*((enf-eof) - de_);
        de_ = enf - eof;
      }
    }

    // if accepted, update
    if (criteria_->accept(lnpMet_, pair_->peTot() + de_, trialType_.c_str(),
                          reject_) == 1) {
      space_->wrap(mpartall_);
      if (!space_->sphereSymMol()) {
        space_->qMolInit(space_->mol()[mpart_.front()]);
      }
      if (dualCut_ > 0) {
        pair_->update(mpart_, 0, "update");
      } else {
        pair_->update(de_);
      }
      if (verbose_ == 1) cout << "accepted " << de_ << std::endl;
      trialAccept();
    } else {
      // skip these rejection steps if nMol == 0
      if (space_->nMol() != 0) {
        space_->restore(mpart_);
        if (space_->cellType() > 0) {
          space_->updateCellofiMol(space_->mol()[mpart_.front()]);
        }
      }
      if (verbose_ == 1) cout << "rejected " << de_ << std::endl;
      trialReject();
    }
  } else if (growType.compare("gcinsdel") == 0) {
    // flip a coin to decide whether to add or delete
    if (uniformRanNum() < 0.5) {
      // add
      trialType_.assign("add");
      mout_("verbose", std::ostringstream().flush() << "attempting to "
            << trialType_);

      // add molecule
      space_->addMol(molType_.c_str());
      pair_->addPart();

      // obtain vector of particle IDs of inserted molecule
      mpart_ = mpartall_ = space_->lastMolIDVec();

      // obtain energy and rosenbluth factors
      double lnwn, en;
      buildAllStep(en, lnwn, 1);
      space_->wrap(mpart_);

      // acceptance criteria
      if (reject_ == 1) {
        de_ = 0;
        lnpMet_ = std::numeric_limits<double>::min();
      } else {
        de_ = en;
        const int iMolIndex = space_->findAddMolListIndex(molType_);
        lnpMet_ = lnwn + log(criteria_->activ(iMolIndex) * space_->vol()
                             / space_->nMol());

        // if dual cut config bias, compute full energy term
        if (dualCut_ > 0) {
          double enf = pair_->multiPartEner(mpart_, 3);
          pair_->update(mpart_, 3, "store");
          lnpMet_ -= criteria_->beta()*(enf - en);
          de_ = enf;
        }
      }

      // if accepted, update
      if (criteria_->accept(lnpMet_, pair_->peTot() + de_, trialType_.c_str(),
                            reject_) == 1) {
        if (dualCut_> 0) {
          pair_->update(mpart_, 3, "update");
        } else {
          pair_->update(de_);
        }
        if (!space_->sphereSymMol()) {
          space_->qMolInit(space_->mol()[mpart_.front()]);
        }
        if (verbose_ == 1) cout << "accepted " << de_ << endl;
        trialAccept();
      } else {
        pair_->delPart(mpart_);
        space_->delPart(mpart_);
        if (verbose_ == 1) cout << "rejected " << de_ << std::endl;
        trialReject();
      }

    // otherwise, delete molecule
    } else {
      trialType_.assign("del");
      mout_("verbose", std::ostringstream().flush() << "attempting to "
            << trialType_);

      double lnwo, eo;
      if (space_->nMol() == 0) {
        reject_ = 1;
      } else {
        // select random molecule to delete
        mpart_ = mpartall_ = space_->randMol();

        // obtain energy and rosenbluth factors
        buildAllStep(eo, lnwo, 0);
      }

      // acceptance criteria
      if (reject_ == 1) {
        de_ = 0;
        lnpMet_ = std::numeric_limits<double>::min();
      } else {
        de_ = -eo;
        const int iMolIndex = space_->findAddMolListIndex(molType_);
        lnpMet_ = -lnwo + log(space_->nMol() / criteria_->activ(iMolIndex)
                              / space_->vol());

        // if dual cut config bias, compute full energy term
        if (dualCut_ > 0) {
          const double eof = pair_->multiPartEner(mpart_, 2);
          pair_->update(mpart_, 2, "store");
          lnpMet_ -= criteria_->beta()*(-eof + eo);
          de_ = -eof;
        }
      }

      // if accepted, update
      if (criteria_->accept(lnpMet_, pair_->peTot() + de_, trialType_.c_str(),
                            reject_) == 1) {
        pair_->delPart(mpart_);
        space_->delPart(mpart_);
        if (dualCut_ > 0) {
          pair_->update(mpart_, 2, "update");
        } else {
          pair_->update(de_);
        }
        if (verbose_ == 1) cout << "accepted " << de_ << std::endl;
        trialAccept();
      } else {
        if (verbose_ == 1) cout << "rejected " << de_ << std::endl;
        trialReject();
      }
    }
  } else {
    ASSERT(0, "unrecognized growType(" << growType << ")");
  }
}

/**
 * implement all steps in cbmc, return energy and rosenbluth factor
 *  newConf == 1, new configuration to compute rosenbluth
 *  newConf != 1, old configuration
 */
void TrialConfigBias::buildAllStep(double &pe, double &lnw, const int newConf) {
  lnw = log(1);
  pe = 0;

  // if avb is on as the first step, update variables
  // tmpart_, mpart_ and region_
  if (stepName_.front().compare("avb") == 0) {
    // check that neighbor list is on a per-atom basis
    ASSERT(pair_->atomCut() != 0,
      "avb in configBias assumes per atom neighbor list");

    // attempt AVB insertion/deletion if gcinsdel
    if (growType.compare("gcinsdel") == 0) {
      lnw = buildAVB(newConf);

    // attempt AVB 2 or 3 if regrowing molecule
    } else if (growType.compare("regrow") == 0) {
      // Because old configurations are computed first, only select
      // targets/regions/etc at that stage
      // For old configurations of CB algorithm with multiple first beads,
      // the opposite region should be used
      if (newConf != 1) {
        // randomly select AVB2 or AVB3
        if (uniformRanNum() < 0.5) {
          lnw = buildAVB2();
        } else {
          lnw = buildAVB3();
        }
      }
    }
  }

  if (reject_ != 1) {
    // label all particles in mpart_ as non physical
    //  each time a particle is placed, label it as physical
    //  at the end, label all particles as physical
    //  this allows intramolecular interactions between particles
    //  assumes mpart was parsed to involve only atoms to be placed
    for (unsigned int im = 0; im < mpart_.size(); ++im) {
      pair_->ipartNotPhysical(mpart_[im]);
    }

    int nAdvSteps = -1e8;
    for (int step = 0; step < nSteps(); step+=nAdvSteps) {
      nAdvSteps = stepOrder_[step];
//      cout << "steptt " << step << " nAdvSteps " << nAdvSteps << endl;
      const int k = stepTrials_[step];
      vector<double> w(k);
      vector<double> en(k);
      vector<double> cpdf(k);
      for (int iTrials = 0; iTrials < k; ++iTrials) {
//        // this line creats a movie
//        pair_->printxyz("asdf", 0);

        // select atoms involved in step, listed in plist, and move it,
        // unless first move in old conf.
        int flag = 1;
        if (newConf != 1) {
          if (iTrials == 0) {
            flag = 0;
          } else {
            flag = 2;
          }
        }

        // place the atoms, recording these in plist
        //  if nAdvSteps>1, do multiple steps simultaneously
        vector<int> plist = buildOneStep(step, flag);
        for (int astep = 1; astep < nAdvSteps; ++astep) {
          vector<int> aplist = buildOneStep(step+astep, flag);
          plist.insert(plist.end(), aplist.begin(), aplist.end());
        }

        // store trial positions
        if (iTrials == 0) {
          space_->xStoreMulti(mpart_, -1);
        } else {
          space_->xStoreMulti(mpart_, -2);
        }

        // label placed atoms (in plist) as physical for energy calculation
        if (iTrials == 0) {
          for (unsigned int ip = 0; ip < plist.size(); ++ip) {
            pair_->ipartIsPhysical(plist[ip]);
          }
        }

        // compute energy and rosenbluth factors
        if (space_->cellType() != 0) {
          space_->updateCellofiMol(space_->mol()[mpart_[0]]);
        }
        if (dualCut_ == 1) pair_->cheapEnergy(1);
        const double pe = pair_->multiPartEner(plist, 1);
        if (dualCut_ == 1) pair_->cheapEnergy(0);
        en[iTrials] = pe;
        w[iTrials] = exp(-criteria_->beta()*pe);
        if (iTrials != 0) w[iTrials] += w[iTrials - 1];
      }
      const double wr = w.back();
      lnw += log(wr/static_cast<double>(k));

      // reject outright if highly unlikely move
      if (wr == 0) {
        reject_ = 1;
        space_->xStoreMulti(mpart_, 0);
      } else {
        int chosen = -1;
        if (newConf == 1) {
          for (int i = 0; i < k; ++i) cpdf[i] = w[i]/wr;
          chosen = ranFromCPDF(cpdf);
        } else {
          chosen = 0;
        }
        space_->xStoreMulti(mpart_, chosen);
        pe += en[chosen];
        // cout << "chosen " << chosen << " " << en[chosen] << endl;
      }
    }
    if (space_->cellType() != 0) {
      space_->updateCellofiMol(space_->mol()[mpart_[0]]);
    }

//    // at the end, check that all particles in mpart_ are labelled as physical
//    for (unsigned int im = 0; im < mpart_.size(); ++im) {
//      ASSERT(pair_->nonphys()[mpart_[im]] == 0,
//        "all particles that are placed must be physical at end of config bias");
//    }
  }
}

/**
 * build the atom in the given step, and return its intermolecular interaction energy
 *  flag == 1, generate new configurations
 *  flag == 2, generate old configurations
 *  flag == 0, return plist, but do not update configuration
 */
vector<int> TrialConfigBias::buildOneStep(const int step, const int flag) {
//  cout << "step name " << stepName_[step] << endl;
  const int iMol = space_->mol()[mpart_[0]];
  const int firstAtom = space_->mol2part()[iMol];
  const int iAtom = stepAtom_[step] + firstAtom;
  vector<int> plist;
  plist.push_back(iAtom);
  if (stepName_[step].compare("freeCOM") == 0) {
    // place iAtom at COM of all other atoms in mpart, for free
    if (flag != 0) space_->setAtomAsCOM(iAtom, mpartall_);
  } else if (stepName_[step].compare("box") == 0) {
    if (flag != 0) {
      // randomly place iAtom in box
      space_->randDisp(iAtom, -1);
    }
  } else if (stepName_[step].compare("avb") == 0) {
    if (flag != 0) {
      if (growType.compare("gcinsdel") == 0) {
        // place iAtom in aggregation volume of tmpart
        space_->avb(iAtom, tmpart_.front(), rAbove_, rBelow_, "in");
      } else if (growType.compare("regrow") == 0) {
        if (flag == 1) {
          // new configuration
          space_->avb(iAtom, tmpart_.front(), rAbove_,
                      rBelow_, region_.c_str());
        } else if (flag == 2) {
          // old configuration
          // string antiRegion("in");
          // if (region_.compare("in") == 0) antiRegion.assign("out");
          // space_->avb(iAtom, tmpart_.front(), rAbove_,
          //  rBelow_, antiRegion.c_str());
          space_->avb(iAtom, tmpartOld_.front(), rAbove_,
                      rBelow_, regionOld_.c_str());
        } else {
          ASSERT(0, "unrecognized flag(" << flag << ") in buildOneStep");
        }
      } else {
        ASSERT(0, "unrecognized growType(" << growType << ") in avb stepping");
      }
    }
  } else if (stepName_[step].compare("bond") == 0) {
    if (flag != 0) {
      // randomly place iAtom in along bond with jAtom
      const int jAtom = stepList_[step][0] + firstAtom;
      double length = stepParam_[step][1];
      const double k = stepParam_[step][0];
      if (fabs(k+1) > doubleTolerance) {
        length = ranBond(criteria_->beta()*k, length, bondExponent);
      }
      space_->setAtomInSphere(iAtom, jAtom, length);
    }
  } else if (stepName_[step].compare("angle") == 0) {
    if (flag != 0) {
      // randomly place iAtom in circle of angle theta w.r.t <ijk
      // and bond length l
      const int jAtom = stepList_[step][0] + firstAtom;
      const int kAtom = stepList_[step][1] + firstAtom;
      double length = stepParam_[step][1];
      double k = stepParam_[step][0];
      if (fabs(k+1) > doubleTolerance) {
        length = ranBond(criteria_->beta()*k, length, bondExponent);
      }
      double theta = stepParam_[step][3];
      k = stepParam_[step][2];

      // rigid angle if k==-1
      if (fabs(k+1) > doubleTolerance) {
        theta = ranAngle(criteria_->beta()*k, theta, angleExponent);
      }
      space_->setAtomInCircle(iAtom, jAtom, kAtom, length, theta);
    }
  } else if (stepName_[step].compare("branch") == 0) {
    // randomly place iAtom and jAtom in branch wrt k and l
    const int jAtom = stepList_[step][0] + firstAtom;
    plist.push_back(jAtom);
    if (flag != 0) {
      const int kAtom = stepList_[step][1] + firstAtom;
      const int lAtom = stepList_[step][2] + firstAtom;
      double Ljk = stepParam_[step][1];
      double Lik = stepParam_[step][3];
      double k = stepParam_[step][0];
      if (fabs(k+1) > doubleTolerance) {
        Ljk = ranBond(criteria_->beta()*k, Ljk, bondExponent);
      }
      k = stepParam_[step][2];
      if (fabs(k+1) > doubleTolerance) {
        Lik = ranBond(criteria_->beta()*k, Lik, bondExponent);
      }
      double Tikl = stepParam_[step][5],
             Tjkl = stepParam_[step][7],
             Tikj = stepParam_[step][9];
      double kikl = stepParam_[step][4], kjkl = stepParam_[step][6],
             kikj = stepParam_[step][8];

      // rigid angles if k==-1
      if ( (fabs(kikl+1) > doubleTolerance) ||
           (fabs(kjkl+1) > doubleTolerance) ||
           (fabs(kikj+1) > doubleTolerance) ) {
        // no branching implemented with mixed flexible/rigid angles
        ASSERT((fabs(kikl + 1) > doubleTolerance) &&
               (fabs(kjkl + 1) > doubleTolerance) &&
               (fabs(kikj + 1) > doubleTolerance),
          "no branching implemented with mixed flexible/rigid angles");

        const double Tikl0 = Tikl, Tjkl0 = Tjkl, Tikj0 = Tikj;

        // alternative - place two vectors in unit sphere,
        // accept by boltzmann weight
        const int maxAttempts = 1e6;
        int accept = 0, i = 0;
        while ( (accept == 0) && (i < maxAttempts) ) {
          vector<double> x1 = ranUnitSphere(3);
          vector<double> x2 = ranUnitSphere(3);
          Tikl = acos(x1[0]);
          Tjkl = acos(x2[0]);
          Tikj = acos(x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2]);
          const double dtikl = Tikl-Tikl0, dtjkl = Tjkl-Tjkl0,
                       dtikj = Tikj-Tikj0;
          if (uniformRanNum() < exp(-criteria_->beta() *
              (kikl*pow(dtikl, angleExponent) + 
               kjkl*pow(dtjkl, angleExponent) +
               kikj*pow(dtikj, angleExponent)))) {
            accept = 1;
          }
          if ( (fabs(kikl) < doubleTolerance) && (dtikl < 0) ) accept = 0;
          if ( (fabs(kjkl) < doubleTolerance) && (dtjkl < 0) ) accept = 0;
          if ( (fabs(kikj) < doubleTolerance) && (dtikj < 0) ) accept = 0;
          ++i;
        }
        ASSERT(i != maxAttempts,
          "maximum attempts reached in selecting branch angles(" << i << ")");

//        // alternative - numerically select values for the angles
//        // do not use this method, it is 'efficient' but does not 
//        // produce the correct distribution
//        double sumAng = 3*PI;
//        const int nMax = 1e4;
//        int i = 0;
//        while (sumAng > 2*PI) {
//          double k = stepParam_[step][4];
//          // rigid angle if k==-1
//          if (fabs(k+1) > doubleTolerance) {
//            Tikl = ranAngle(criteria_->beta()*k, Tikl0, angleExponent);
//          }
//          k = stepParam_[step][6];
//          // rigid angle if k==-1
//          if (fabs(k+1) > doubleTolerance) {
//            Tjkl = ranAngle(criteria_->beta()*k, Tjkl0, angleExponent);
//          }
//          k = stepParam_[step][8];
//          // rigid angle if k==-1
//          if (fabs(k+1) > doubleTolerance) {
//            Tikj = ranAngle(criteria_->beta()*k, Tikj0, angleExponent);
//          }
//          sumAng = Tikl + Tjkl + Tikj;
//          ++i;
//
//          // also, check that one angle isn't greater than
//          // the other two combined.
//          if ( (Tikl/sumAng > 0.5) ||
//               (Tjkl/sumAng > 0.5) || (Tikj/sumAng > 0.5) ) sumAng = 3*PI;
//          ASSERT(i <= nMax, "reached nMax in finding a set of 3 angles");
//        }
        
        // end of two alternatives

        // cout << "Tikl " << cos(Tikl) << " " << Tikl/PI*180 << " Tjkl "
        //  << cos(Tjkl) << " " << Tjkl/PI*180 << " Tikj " << cos(Tikj)
        //  << " " << Tikj/PI*180 << " sum " << sumAng << " "
        //  << sumAng/PI*180 << endl;
      }
      space_->setAtomInCircle(jAtom, kAtom, lAtom, Ljk, Tjkl);
      space_->setAtomInBranch(lAtom, jAtom, iAtom, kAtom, Tikl, Tikj, Lik);
    }
  } else {
    ASSERT(0, "unrecognized stepName(" << stepName_[step] << ") in molType("
           << molType_ << ") called in stepper");
  }
  return plist;
}

/*
 * reduce list of particles for a molecule to only particles involved
 * in CB move
 */
void TrialConfigBias::mpartCBparse() {
  vector<int> plist;
  for (int step = 0; step < nSteps(); ++step) {
    vector<int> plistStep = buildOneStep(step, 0);
    plist.insert(plist.end(), plistStep.begin(), plistStep.end());
  }
  std::sort(plist.begin(), plist.end());
  mpart_ = plist;
}

/**
 * For AVB, convert tmpart from list of all sites in a target molecule
 * to a random site available for AVB. Also return neighbors of that same type
 */
vector<int> TrialConfigBias::mpart2neigh(vector<int> &tmpart) {
  // list atoms in tmpart that are available for avb, provided in stepList
  const int tMol = space_->mol()[tmpart.front()];
  tmpart = stepList_.front();
  for (unsigned int i = 0; i < tmpart.size(); ++i) {
    tmpart[i] += space_->mol2part()[tMol];
  }

  // select a random atom in tmpart
  const int index = uniformRanNum(0, static_cast<int>(tmpart.size())-1);
  const int jpart = tmpart[index];
  tmpart.clear();
  tmpart.push_back(jpart);

  // list neighbors of jpart that are of the same type of atom as the stepAtom
  vector<int> tneigh = pair_->neigh()[jpart];
  for (int i = static_cast<int>(tneigh.size()) - 1; i >= 0; --i) {
    const int jMol = space_->mol()[tneigh[i]];
    if (tneigh[i] != (space_->mol2part()[jMol] + stepAtom_.front())) {
      tneigh.erase(tneigh.begin() + i);
    }
  }
  return tneigh;
}

/**
 * initialize buildAllStep for AVB insertions/deletions, return prefactor
 */
double TrialConfigBias::buildAVB(const int newConf) {
  double lnw = log(1);
  // reject if less than two molecules present
  if (space_->nMol() <= 1) {
    reject_ = 1;
  } else {
    if (newConf == 1) {
      // if adding, select target of different molecule
      const int iMol = space_->mol()[mpart_.front()];
      tmpart_ = space_->randMolDiff(iMol);
    } else {
      // if deleting, select random molecule
      tmpart_ = space_->randMol();
    }

    // convert tmpart to list target atom, and return neighbors of that atom
    vector<int> tneigh = mpart2neigh(tmpart_);
    int nIn = static_cast<int>(tneigh.size());
    if (newConf == 1) ++nIn;
    if (nIn == 0) {
      reject_ = 1;   // reject deletion move because no particles present
    } else {
      lnw = log(vIn_*static_cast<double>(space_->nMol()-1)
                /static_cast<double>(space_->nMol())/static_cast<double>(nIn));

      // cancel out assumed prefactor for insertion into entire box
      lnw += log(static_cast<double>(space_->nMol()) / space_->vol());

      // if deleting, pick an atom from tneigh, and select that molecule
      if (newConf != 1) {
        const int kpart = tneigh[uniformRanNum(0, tneigh.size()-1)];
        mpart_ = mpartall_ = space_->imol2mpart(space_->mol()[kpart]);
      }
    }
  }
  return lnw;
}

/**
 * initialize buildAllStep for AVB2 moves, return prefactor
 *
 *   AVB2 Algorithm of Chen and Siepmann,
 *     "Improving the Efficiency of the Aggregation Volume Bias Monte Carlo
 *      Algorithm"
 *     J. Phys. Chem. B, 2001, 105 (45), pp 11275-11282
 */
double TrialConfigBias::buildAVB2() {
  double lnw = log(1);

  // reject if less than two molecules present
  if (space_->nMol() <= 1) {
    reject_ = 1;
  } else {
    // For a given regrow, must choose whether to do in->out or out->in
    // according to AVB2
    // This choice will be recorded in region_ (e.g., region_ == in
    // for out->in move)
    // probabily of choosing in->out or out->in is hardcoded at 50%
    const double pBias = 0.5;
    // select a random molecule, tmpart, as target for AVB2 swap move
    tmpart_ = space_->randMol();
    // convert tmpart to target atom, and return neighbors of that atom
    vector<int> tneigh = mpart2neigh(tmpart_);
    tmpartOld_ = tmpart_;
    const int nIn = static_cast<int>(tneigh.size());
    vector<int> tNotNeigh = iAtom2notNeigh(tmpart_.front(), tneigh);
    const int nOut = static_cast<int>(tNotNeigh.size());
    const double vOut = space_->vol() - vIn_;

    // choose in->out or out->in
    if (uniformRanNum() < pBias) {
      // perform an out->in AVB2 regrow move
      region_.assign("in");
      regionOld_.assign("out");

      // randomly select an atom of type stepAtom in the out region of tmpart
      if (nOut == 0) {
        reject_ = 1;
      } else {
        // randomly choose atom in tNotNeigh to determine mpart_
        chooseMPart(tNotNeigh);

        // obtain the prefactor
        lnw = log((1. - pBias) * vIn_ * nOut / pBias / vOut / (nIn + 1));
        lnw *= -1;   // invert because old configuration
      }
    } else {
      // perform an in->out AVB2 regrow move
      region_.assign("out");
      regionOld_.assign("in");
      if (nIn == 0) {
        reject_ = 1;
      } else {
        // pick a random atom from tneigh to define mpart_
        chooseMPart(tneigh);

        // obtain the prefactor
        lnw = log(pBias * vOut * nIn / (1. - pBias) / vIn_ / (nOut + 1));
        lnw *= -1;   // invert because old configuration
      }
    }
  }
  return lnw;
}

/**
 * initialize buildAllStep for AVB3 moves, return prefactor
 *
 *   AVB3 Algorithm of Chen and Siepmann,
 *     "Improving the Efficiency of the Aggregation Volume Bias Monte Carlo
 *      Algorithm"
 *     J. Phys. Chem. B, 2001, 105 (45), pp 11275-11282
 */
double TrialConfigBias::buildAVB3() {
  double lnw = log(1);

  // reject if less than three molecules present
  if (space_->nMol() <= 2) {
    reject_ = 1;
  } else {
    // For a given regrow, must choose whether to do in->out or out->in
    //  This choice will be recorded in region_ (e.g., region_ == in
    //  for out->in move)
    //  probabily of choosing in-->out or out-> is hardcoded at 50%
    const double pBias = 0.5;

    // select random molecules j and k(!=j) with non-overlapping bonded
    // regions (rij>2rabove) as targets for swap move
    vector<int> kmpart, jmpart = space_->randMol();
    const vector<int> jneigh = mpart2neigh(jmpart);
    const int nInJ = static_cast<int>(jneigh.size());

    // select kmpart != jmpart
    vector<int> kneigh;
    int bonded = 1, randPicks = 0;
    const int randPicksMax = 2 * space_->nMol();
    while ( (bonded == 1) && (randPicks != randPicksMax) ) {
      kmpart = space_->randMolDiff(space_->mol()[jmpart.front()]);
      kneigh = mpart2neigh(kmpart);
      bonded = space_->bonded(jmpart[0], kmpart[0], 2.*rAbove_, 0.);
      ++randPicks;
    }

    // avb3 fails if it cannot find j and k with non-overlapping bonded
    // regions
    if (randPicks == randPicksMax) {
      reject_ = 1;
    } else {
      const int nInK = static_cast<int>(kneigh.size());
      const double vOut = space_->vol() - vIn_;

      // select move with probability pbias
      if (uniformRanNum() < pBias) {
        // select one particle either from in region of k or out region
        // of j (equal prob) to go to in region of j
        if (uniformRanNum() < 0.5) {
          // select particle from in region of k
          if (nInK == 0) {
            reject_ = 1;
          } else {
            chooseMPart(kneigh);
            // in k -> in j
            lnw = log((1 - pBias)*vIn_*nInK/pBias/vIn_/(nInJ + 1));
            lnw *= -1;   // invert because old configuration
            regionOld_.assign("in");
            tmpartOld_ = kmpart;
            if (verbose_ == 1) {
              cout << "selected particle " << mpart_[0]
                   << " from in region of k " << kmpart[0] << endl;
            }
          }
        } else {
          // select particle from out region of j
          vector<int> jNotNeigh = iAtom2notNeigh(jmpart.front(), jneigh);
          const int nOutJ = static_cast<int>(jNotNeigh.size());
          if (nOutJ == 0) {
            reject_ = 1;
          } else {
            chooseMPart(jNotNeigh);
            // out j -> in j
            lnw = log((1 - pBias)*vIn_*nOutJ/pBias/vOut/(nInJ + 1));
            lnw *= -1;   // invert because old configuration
            regionOld_.assign("out");
            tmpartOld_ = jmpart;
          }
        }

        // move mpart to the in region of j
        region_.assign("in");
        tmpart_ = jmpart;
      } else {
        // select mpart particle in the in region of j
        if (nInJ == 0) {
          reject_ = 1;
        } else {
          chooseMPart(jneigh);
          regionOld_.assign("in");
          tmpartOld_ = jmpart;

          // move mpart either to in region of k or to out region
          // of j (equal prob)
          if (uniformRanNum() < 0.5) {
            // move mpart to the in region of k
            region_.assign("in");
            tmpart_ = kmpart;
            // in j -> in k
            lnw = log(pBias*vIn_*nInJ/(1 - pBias)/vIn_/(nInK + 1));
            lnw *= -1;   // invert because old configuration

          } else {
            // move mpart to the out region of j
            region_.assign("out");
            tmpart_ = jmpart;
            vector<int> jNotNeigh = iAtom2notNeigh(jmpart.front(), jneigh);
            const int nOutJ = static_cast<int>(jNotNeigh.size());
            // in j -> out j
            lnw = log(pBias*vOut*nInJ/(1 - pBias)/vIn_/(nOutJ + 1));
            lnw *= -1;   // invert because old configuration
          }
        }
      }
    }
  }
  return lnw;
}

/**
 * randomly choose an atom in list to define mpart
 */
void TrialConfigBias::chooseMPart(const vector<int> &list) {
  const int index = uniformRanNum(0, static_cast<int>(list.size())-1);
  mpart_ = space_->imol2mpart(space_->mol()[list[index]]);
  mpartall_ = mpart_;
}

/**
 * return list of atoms in the out region mpart
 */
vector<int> TrialConfigBias::iAtom2notNeigh(
  const int iAtom, const vector<int> &neigh) {
  // list all neighbor atoms of same type as iAtom, including iAtom
  vector<int> neighSelf = neigh;
  neighSelf.push_back(iAtom);
  const int type = space_->type()[iAtom];

  // generate list of all atoms of same type as iAtom that are not neighbors
  vector<int> notNeigh;
  for (int jAtom = 0; jAtom < space_->natom(); ++jAtom) {
    if (space_->type()[jAtom] == type) {
      if (!findInList(jAtom, neighSelf)) {
        notNeigh.push_back(jAtom);
      }
    }
  }
  return notNeigh;
}

/**
 * initialize trials in given MC class based on reading input file
 */
void TrialConfigBias::file2MC(const char* fileName, MC* mc) {
  // use file extension to determine whether or not to use JSON or LMP
  // data files
  if (myTrim(".", fileName) == "json") {
    #ifdef JSON_
      json2MC(fileName, mc);
    #else
      ASSERT(0, "JSON files are not implemented, but file " << fileName <<
        " was attempted to be read");
    #endif  // JSON_

  } else {
    lmp2MC(fileName, mc);
  }
}

/**
 * initalize Configurational Bias moves, given above flags, from LMP data file
 */
void TrialConfigBias::lmp2MC(const char* fileName, MC* mc) {
  // open LAMMPS data file
  std::ifstream file(fileName);
  ASSERT(file.good(), "cannot find lammps DATA file " << fileName);

  // read until Configuration section
  readUntil("ConfigurationalBias", file);

  // read number of types of moves
  int cbmctypes;
  file >> cbmctypes;

  // read each type
  for (int itype = 0; itype < cbmctypes; ++itype) {
    string descript, line;
    getline(file, line);
    int nsteps;
    double weight;
    file >> descript >> nsteps >> weight;
    mc->weight = weight;
    shared_ptr<TrialConfigBias> trial = make_shared<TrialConfigBias>(*this);

    // skip gcinsdel moves if not using grand canonical ensemble
    bool skip = false;
    if ( (insDelFlag == 0) && (descript.compare("gcinsdel") == 0) ) skip = true;
    if (insDelFlag == 2) descript.assign("regrow");
    if (!skip) {
      trial->growType = descript;
      mc->initTrial(trial);
    }
    for (int istep = 0; istep < nsteps; ++istep) {
      getline(file, line);
      int iAtom, ntrials;
      file >> iAtom >> ntrials >> descript;
      iAtom -= 1;
      if (descript.compare("bond") == 0) {
        int jAtom;
        file >> jAtom;
        jAtom -= 1;
        if (!skip) trial->addStepBond(iAtom, jAtom, ntrials);
      } else if (descript.compare("angle") == 0) {
        int jAtom, kAtom;
        file >> jAtom >> kAtom;
        jAtom -= 1; kAtom -= 1;
        if (!skip) trial->addStepAngle(iAtom, jAtom, kAtom, ntrials);
      } else if (descript.compare("branch") == 0) {
        int jAtom, kAtom, lAtom;
        file >> jAtom >> kAtom >> lAtom;
        jAtom -= 1; kAtom -= 1; lAtom -= 1;
        if (!skip) trial->addStepBranch(iAtom, jAtom, kAtom, lAtom, ntrials);
      } else if (descript.compare("box") == 0) {
        if (!skip) trial->addStepBox(iAtom, ntrials);
      } else if (descript.compare("freeCOM") == 0) {
        if (!skip) trial->addStepFreeCOM(iAtom);
      } else if (descript.compare("avb") == 0) {
        double rAbove, rBelow;
        int natom;
        vector<int> tmpart;
        file >> rAbove >> rBelow >> natom;
        for (int ipart = 0; ipart < natom; ++ipart) {
          int part;
          file >> part;
          tmpart.push_back(part-1);
        }
        if (!skip) {
          trial->addStepAVB(iAtom, tmpart, rAbove, rBelow, ntrials);
          mc->neighAVBInit(trial->rAbove(), trial->rBelow());
        }
      } else {
        ASSERT(0, "unrecognized configurational bias move type(" << descript
          << ") when reading file(" << fileName << ")");
      }
    }
    if (!skip) {
      trial->updateStepParams();
      trial->initDualCut(dualCut_);
    }
  }
}


/**
 * initalize Configurational Bias moves, given above flags, from JSON
 * data file
 */
void TrialConfigBias::json2MC(const char* fileName, MC* mc) {
#ifdef JSON_
  // open a JSON file
  std::ifstream file(fileName);
  ASSERT(file.good(), "cannot find json DATA file " << fileName);
  json j;
  file >> j;

  // for each configurational bias move set
  for (unsigned int icb = 0; icb < j["configBias"].size(); ++icb) {
    json jset = j["configBias"][icb];
    string descript = jset["type"];
    double weight;
    weight = jset["weight"];
    mc->weight = weight;
    shared_ptr<TrialConfigBias> trial = make_shared<TrialConfigBias>(*this);

    // skip gcinsdel moves if not using grand canonical ensemble
    bool skip = false;
    if ( (insDelFlag == 0) && (descript.compare("gcinsdel") == 0) ) skip = true;
    if (insDelFlag == 2) descript.assign("regrow");
    if (!skip) {
      trial->growType = descript;
      mc->initTrial(trial);
    }

    // for each step in the move set
    for (unsigned int istep = 0; istep < jset["steps"].size(); ++istep) {
      json jstep = jset["steps"][istep];
      descript = jstep["type"];
      const int ntrials = jstep["trials"];
      int iAtom = jstep["atom"][0];
      iAtom -= 1;
      if (descript.compare("bond") == 0) {
        int jAtom = jstep["bondTo"];
        jAtom -= 1;
        if (!skip) trial->addStepBond(iAtom, jAtom, ntrials);
      } else if (descript.compare("angle") == 0) {
        int jAtom = jstep["bondTo"];
        int kAtom = jstep["angleRef"];
        jAtom -= 1; kAtom -= 1;
        if (!skip) trial->addStepAngle(iAtom, jAtom, kAtom, ntrials);
      } else if (descript.compare("branch") == 0) {
        int jAtom = jstep["atom"][1];
        int kAtom = jstep["bondTo"];
        int lAtom = jstep["angleRef"];
        jAtom -= 1; kAtom -= 1; lAtom -= 1;
        if (!skip) trial->addStepBranch(iAtom, jAtom, kAtom, lAtom, ntrials);
      } else if (descript.compare("box") == 0) {
        if (!skip) trial->addStepBox(iAtom, ntrials);
      } else if (descript.compare("freeCOM") == 0) {
        if (!skip) trial->addStepFreeCOM(iAtom);
      } else if (descript.compare("avb") == 0) {
        ASSERT(0, "avb not implemented in JSON");
//        double rAbove, rBelow;
//        int natom;
//        vector<int> tmpart;
//        file >> rAbove >> rBelow >> natom;
//        for (int ipart = 0; ipart < natom; ++ipart) {
//          int part;
//          file >> part;
//          tmpart.push_back(part-1);
//        }
//      if (!skip) lastCB().addStepAVB(iAtom, tmpart, rAbove, rBelow, ntrials);
//        mc->neighAVBInit(mc->lastCB().rAbove(), mc->lastCB().rBelow());
      } else {
        ASSERT(0, "unrecognized configurational bias move type(" << descript
          << ") when reading file(" << fileName << ")");
      }
    }
    if (!skip) {
      trial->updateStepParams();
      trial->initDualCut(dualCut_);
    }
  }
#endif  // JSON_
}

/**
 * initialize dual-cut configurational bias
 */
void TrialConfigBias::initDualCut(const int flag) {
  dualCut_ = flag;
  if (dualCut_ == 1) {
    // set cell list according to largest sigma
    space_->updateCells(pair_->sigMax());
  } else if (dualCut_ == 0 || dualCut_ == 2) {
    // do nothing
  } else {
    ASSERT(0, "unrecognized flag for dualCut(" << dualCut_ << ")");
  }
}

}  // namespace feasst

