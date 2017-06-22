/**
 * \file
 *
 * \brief pairwise interactions
 *
 * Implementation of the pair class
 */

#include "./pair.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

Pair::Pair(Space* space,
  const double rCut)   //!< interaction cut-off distance
  : space_(space),
    rCut_(rCut) {
  className_.assign("Pair");
  defaultConstruction();
}
Pair::Pair(Space* space, const char* fileName)
  : space_(space) {
  className_.assign("Pair");
  ASSERT(fileExists(fileName), "restart file(" << fileName
         << ") doesn't exist");
  rCut_ = fstod("rCut", fileName);
  defaultConstruction();

  // cout << "initialize random number generator" << endl;
  initRNG(fileName);

  // cout << "initialize long range corrections" << endl;
  string strtmp = fstos("lrcFlag", fileName);
  if (!strtmp.empty()) {
    lrcFlag = atoi(strtmp.c_str());
  }

  // cout << "initialize data files" << endl;
  strtmp = fstos("nDataFiles", fileName);
  if (!strtmp.empty()) {
    const int ndf = stoi(strtmp);
    for (int i = 0; i < ndf; ++i) {
      stringstream ss;
      ss << "DATAFILE" << i;
      const string datastr = fstos(ss.str().c_str(), fileName);
      // cout << "df " << i << " file " << lmpstr << endl;
      initData(datastr.c_str());
    }
  }

  // cout << "initialize order" << endl;
  strtmp = fstos("myuniqueorderParam", fileName);
  if (!strtmp.empty()) {
    orderOn_ = 1;
    order_ = stod(strtmp);
    orderMax_ = fstod("myuniqueorderMax", fileName);
    orderMin_ = fstod("myuniqueorderMin", fileName);
    strtmp = fstos("myuniqueorderName", fileName);
    if (!strtmp.empty()) orderName_ = strtmp;
    setOrder(order_);
  }

  // cout << "initialize manually set epsij" << endl;
  strtmp = fstos("nepsijsetRecord", fileName);
  if (!strtmp.empty()) {
    const int nes = stoi(strtmp);
    for (int i = 0; i < nes; ++i) {
      vector<double> eps;
      for (int j = 0; j < 3; ++j) {
        stringstream ss;
        ss << "epsijsetRecord" << i << j;
        eps.push_back(fstod(ss.str().c_str(), fileName));
      }
      epsijset(feasstRound(eps[0]), feasstRound(eps[1]), eps[2]);
    }
  }

  // cout << "set rcutij" << endl;
  rCutij_.resize(epsij_.size(), vector<double>(epsij_.size(), rCut_));
  for (unsigned int i = 0; i < epsij_.size(); ++i) {
    for (unsigned int j = i; j < epsij_.size(); ++j) {
      stringstream ss;
      ss << "rCut" << i << j;
      string strtmp = fstos(ss.str().c_str(), fileName);
      if (!strtmp.empty()) {
        const double rc = stod(strtmp);
        rCutijset(i, j, rc);
      }
    }
  }

  // cout << "initialize neighbor variables" << endl;
  strtmp = fstos("atomCut", fileName);
  if (!strtmp.empty()) {
    atomCut_ = atoi(strtmp.c_str());
    initAtomCut(atomCut_);
  }
  strtmp = fstos("neighTypeScreen", fileName);
  if (!strtmp.empty()) {
    neighTypeScreen_ = atoi(strtmp.c_str());
    if (neighTypeScreen_ == 1) {
      neighTypeSet(fstoi("neighTypeIndex0", fileName));
      neighTypeSet(fstoi("neighTypeIndex1", fileName));
    }
  }

  strtmp = fstos("intraMolecInteract", fileName);
  if (!strtmp.empty()) {
    intra_ = atoi(strtmp.c_str());
  }

  strtmp = fstos("sigrefFlag", fileName);
  if (!strtmp.empty()) {
    sigrefFlag_ = atoi(strtmp.c_str());
  }
}

/**
 * defaults in constructor
 */
void Pair::defaultConstruction() {
  verbose_ = 0;
  orderOn_ = 0;
  order_ = -1;
  orderMin_ = order_;
  orderMax_ = order_;
  rCutSq_ = rCut_ * rCut_;
  fill(rCut_, rCutij_);
  rCutMaxAll_ = rCut_;
  dimen_ = space_->dimen(),
  f_.resize(space_->natom(), vector<double>(dimen_));
  vr_.resize(space_->natom(), vector<vector<double> >(
    dimen_, vector<double>(dimen_)));
  pe_.resize(space_->natom());
  lrcFlag = 1;    // compute long range corrections by default
  sigepsDefault_ = 1;
  eps_.resize(1); eps_.at(0) = 1.;
  sig_.resize(1); sig_.at(0) = 1.;
  neighOn_ = false;
  atomCut_ = 1;
  fastDel_ = false;
  fastDelMol_ = -1;
  ASSERT((space_->minl() == 0) || (rCut_ <= space_->minl()/2.),
    "the selected pair-wise interaction cutoff, rCut(" << rCut_
    << ") is too large for the minimum box length(" << space_->minl()/2.
    << ").");
  cheapEnergy_ = false;
  if (atomCut_ == 0) {
    neigh_.resize(space_->nMol(), vector<int>(0));
    neighCut_.resize(space_->nMol(), vector<int>(0));
  } else {
    neigh_.resize(space_->natom(), vector<int>(0));
    neighCut_.resize(space_->natom(), vector<int>(0));
  }
  neighTypeScreen_ = 0;
  neighCutOn_ = 0;
  peMapOn_ = 0;
  intra_ = 0;
  addPart();
  sigrefFlag_ = 0;
}

/**
 * reset object pointers
 */
void Pair::reconstruct(Space* space) {
  space_ = space;
  Base::reconstruct();
}

/**
 * write restart file
 */
void Pair::writeRestartBase(const char* fileName) {
  fileBackUp(fileName);
  std::ofstream file(fileName);
  file << "# className " << className_ << endl;
  if (lrcFlag != 1) file << "# lrcFlag " << lrcFlag << endl;
  if (orderOn_ == 1) {
    file << "# myuniqueorderParam " << order_ << endl;
    file << "# myuniqueorderMax " << orderMax_ << endl;
    file << "# myuniqueorderMin " << orderMin_ << endl;
    if (!orderName_.empty()) {
      file << "# myuniqueorderName " << orderName_ << endl;
    }
  }
  file << std::setprecision(std::numeric_limits<double>::digits10+2)
       << "# rCut " << rCut_ << endl;
  if (initDataFiles_.size() != 0) {
    file << "# nDataFiles " << initDataFiles_.size() << endl;
    for (unsigned int i = 0; i < initDataFiles_.size(); ++i) {
      file << "# DATAFILE" << i << " " << initDataFiles_[i] << endl;
    }
  }
  const int nes = static_cast<int>(epsijsetRecord_.size());
  if (nes != 0) {
    file << "# nepsijsetRecord " << nes << endl;
    for (int i = 0; i < nes; ++i) {
      for (int j = 0; j < 3; ++j) {
        file << std::setprecision(std::numeric_limits<double>::digits10+2)
             << "# epsijsetRecord" << i << j << " " << epsijsetRecord_[i][j]
             << endl;
      }
    }
  }

  // print rCutij
  if (rCutij_.size() > 0) {
    for (unsigned int i = 0; i < epsij_.size(); ++i) {
      for (unsigned int j = i; j < epsij_.size(); ++j) {
        if (rCutij_[i][j] != rCut_) {
          file << std::setprecision(std::numeric_limits<double>::digits10+2)
               << "# rCut" << i << j << " " << rCutij_[i][j] << endl;
        }
      }
    }
  }

  // write neighbor list variables
  file << "# atomCut " << atomCut_ << endl;
  if (neighTypeScreen_ != 0) {
    file << "# neighTypeScreen " << neighTypeScreen_ << endl;
    file << "# neighTypeIndex0 " << neighType_[0] << endl;
    file << "# neighTypeIndex1 " << neighType_[1] << endl;
  }

  if (intra_ != 0) file << "# intraMolecInteract " << intra_ << endl;
  if (sigrefFlag_ != 0) file << "# sigrefFlag " << sigrefFlag_ << endl;

  // write random number generator state
  writeRngRestart(fileName);
}

/**
 * Returns the pressure of the system
 */
double Pair::pressure(const double beta     //!< inverse temperature
  ) {
  return (space_->nMol()/beta + vrTot()/3.)/space_->vol();
}

/**
 * delete particle ipart
 */
void Pair::delPartBase(const int ipart) {  //!< particle number to delete
  // error check that particle exists
  ASSERT(ipart < space_->natom(), "cannot delete particle that does not exist,"
    << "ipart: " << ipart << " when there are only natom: " << space_->natom());

  f_.erase(f_.begin() + ipart);
  pe_.erase(pe_.begin() + ipart);
  vr_.erase(vr_.begin() + ipart);

  // update neighbor list
  if (neighOn_) {
    if ( (atomCut_ == 0) || (space_->natom() == space_->nMol()) ) {
      // check if deleting first atom in molecule
      int iMol = space_->mol()[ipart];
      if (space_->mol2part()[iMol] == ipart) {
        vector<int> mpart(1, iMol);
        eraseNeigh_(mpart, &neigh_);
        eraseNeigh_(mpart, &neighCut_);
      }
    }
  }
}

/**
 * delete particles mpart
 */
void Pair::delPartBase(const vector<int> mpart   //!< particle numbers to delete
  ) {
  fastDel_ = space_->fastDelApplicable(mpart);
  if (fastDel_) fastDelMol_ = space_->mol()[mpart.front()];
  for (int i = static_cast<int>(mpart.size()) - 1; i >= 0; --i) {
    delPartBase(mpart[i]);
  }
  if (neighOn_ && (atomCut_ == 1) && (space_->natom() != space_->nMol())) {
    eraseNeigh_(mpart, &neigh_);
    eraseNeigh_(mpart, &neighCut_);
  }
}

/**
 * erase mpart from neighlist
 *  if atomCut is 1, then iMol and jMol are atom identifiers, not molecules
 */
void Pair::eraseNeigh_(
  const vector<int> mpart,        //!< molecules(or atoms) to erase
  vector<vector<int> > *neighPtr) {  //!< neighborlist to update
  vector<vector<int> >& neigh = *neighPtr;
  if (neigh.size() > 0) {
    if (fastDel_) {
      for (int impart = static_cast<int>(mpart.size()) - 1; impart >= 0;
           --impart) {
        const int iMol = mpart[impart];
        // cout << "erasing neigh in fastdel for iMol " << iMol << endl;
        int jMol = space_->nMol() - 1;
        if (atomCut_ == 1) {
          jMol = space_->natom() - static_cast<int>(mpart.size()) + impart;
        }
        neigh[iMol] = neigh[jMol];
        neigh.erase(neigh.begin() + jMol);
        for (unsigned int i = 0; i < neigh.size(); ++i) {
          // remove iMol
          for (unsigned int j = 0; j < neigh[i].size(); ++j) {
            if (neigh[i][j] == iMol) {
              neigh[i].erase(neigh[i].begin() + j);
            }
          }

          // swap jMol into iMol
          for (unsigned int j = 0; j < neigh[i].size(); ++j) {
            if (neigh[i][j] == jMol) {
              neigh[i][j] = iMol;
            }
          }
        }
      }
    } else {
      for (int impart = mpart.size() - 1; impart >= 0; --impart) {
        const int iMol = mpart[impart];
        neigh.erase(neigh.begin() + iMol);
        for (unsigned int i = 0; i < neigh.size(); ++i) {
          // remove iMol
          for (unsigned int j = 0; j < neigh[i].size(); ++j) {
            if (neigh[i][j] == iMol) {
              neigh[i].erase(neigh[i].begin() + j);
            }
          }

          // count down iMol tags
          if (fastDel_ == false) {
            for (unsigned int j = 0; j < neigh[i].size(); ++j) {
              if (neigh[i][j] >= iMol) {
                --neigh[i][j];
              }
            }
          }
        }
      }
    }
  }
}

/**
 * add particle ipart
 */
void Pair::addPartBase() {
  f_.resize(space_->natom(), vector<double>(space_->dimen()));
  vr_.resize(space_->natom(), vector<vector<double> >(
    space_->dimen(), vector<double>(space_->dimen())));
  pe_.resize(space_->natom());
  nonphys_.resize(space_->natom(), 0);
  if (atomCut_ == 0) {
    neigh_.resize(space_->nMol(), vector<int>(0));
    neighCut_.resize(space_->nMol(), vector<int>(0));
  } else {
    neigh_.resize(space_->natom(), vector<int>(0));
    neighCut_.resize(space_->natom(), vector<int>(0));
  }
}

/**
 * initialize the pair parameters eps and sig
 */
void Pair::initPairParam(
  const vector<double> eps,       //!< epsilon pair parameter
  const vector<double> sig,       //!< sigma pair parameter
  const vector<double> sigref) {  //!< sigma pair parameter
  const int ntype = eps.size();
  ASSERT(sig.size() == eps.size(),
    "input parameters eps and sig in initPairParm must be of equal size");

  sigepsDefault_ = 0;

  eps_.resize(ntype);
  sig_.resize(ntype);
  sigRef_.resize(ntype);
  for (int i = 0; i < ntype; ++i) {
    eps_[i] = eps[i];
    sig_[i] = sig[i];
    if (sigrefFlag_ == 1) sigRef_[i] = sigref[i];
  }

  epsij_.resize(ntype, vector<double>(ntype));
  sigij_.resize(ntype, vector<double>(ntype));
  sigRefij_.resize(ntype, vector<double>(ntype));
  for (int i = 0; i < ntype; ++i) {
    epsij_[i].resize(ntype);
    sigij_[i].resize(ntype);
    if (sigrefFlag_ == 1) sigRefij_[i].resize(ntype);
  }
  for (int i = 0; i < ntype; ++i) {
    for (int j = 0; j < ntype; ++j) {
      epsij_[i][j] = sqrt(eps_[i]*eps_[j]);
      sigij_[i][j] = (sig_[i]+sig_[j])/2;
      if (sigrefFlag_ == 1) sigRefij_[i][j] = (sigRef_[i]+sigRef_[j])/2;
    }
  }
}

/**
 * initialize neighbor list
 */
void Pair::initNeighList(
  const double neighAbove,    //!< upper cut off distance for neighbor list
  const double neighBelow) {  //!< lower cut off distance for neighbor list
  neighAbove_ = neighAbove;
  neighAboveSq_ = neighAbove*neighAbove;
  neighBelow_ = neighBelow;
  neighBelowSq_ = neighBelow*neighBelow;
  neighOn_ = true;
  buildNeighList();
}

/**
 * build full neighbor list of each molecule
 */
void Pair::buildNeighList() {
  // shorthand for read-only space variables
  const int dimen = space_->dimen();
  int nMol = space_->nMol();
  const vector<double> x = space_->x();
  const vector<double> l = space_->l();
  const vector<int> mol2part = space_->mol2part();
  const vector<int> mol = space_->mol();
  const vector<int> type = space_->type();

  // empty neighbor list and resize
  neigh_.clear();
  neighCut_.clear();
  if (atomCut_ == 0) {
    neigh_.resize(space_->nMol(), vector<int>(0));
    neighCut_.resize(space_->nMol(), vector<int>(0));
  } else {
    nMol = space_->natom();
    neigh_.resize(space_->natom(), vector<int>(0));
    neighCut_.resize(space_->natom(), vector<int>(0));
  }

  // loop through pairs of molecules or atoms
  for (int iMol = 0; iMol < nMol - 1; ++iMol) {
    int ipart;
    if (atomCut_ == 0) {
      ipart = mol2part[iMol];
    } else {
      ipart = iMol;
    }
    for (int jMol = iMol + 1; jMol < nMol; ++jMol) {
      if ( (atomCut_ == 0) || (mol[iMol] != mol[jMol]) ) {
        int jpart;
        if (atomCut_ == 0) {
          jpart = mol2part[jMol];
        } else {
          jpart = jMol;
        }

        // separation distance with periodic boundary conditions
        vector<double> rij(space_->dimen());
        for (int dim = 0; dim < dimen; ++dim) {
          rij[dim] = x[dimen*ipart+dim] - x[dimen*jpart+dim];
        }
        const vector<double> dx = space_->pbc(rij);
        for (int dim = 0; dim < dimen; ++dim) {
          rij[dim] += dx[dim];
        }
        const double r2 = vecDotProd(rij, rij);

        if (r2 < rCutSq_) {
          if (neighTypeScreen_ == 1) {
            if ( ( (neighType_[0] == type[ipart]) &&
                   (neighType_[1] == type[jpart]) ) ||
                 ( (neighType_[1] == type[ipart]) &&
                   (neighType_[0] == type[jpart]) ) ) {
              neighCut_[iMol].push_back(jMol);
              neighCut_[jMol].push_back(iMol);
              if ( (r2 < neighAboveSq_) && (r2 > neighBelowSq_) ) {
                neigh_[iMol].push_back(jMol);
                neigh_[jMol].push_back(iMol);
              }
            }
          } else {
            neighCut_[iMol].push_back(jMol);
            neighCut_[jMol].push_back(iMol);
            if ( (r2 < neighAboveSq_) && (r2 > neighBelowSq_) ) {
              neigh_[iMol].push_back(jMol);
              neigh_[jMol].push_back(iMol);
            }
          }
        }
      }
    }
  }
}

/**
 * stores, restores or updates neighbor list variables to avoid recompute of
 * entire configuration after every change
 *  flag==0, energy of old configuration (no moves, ins or dels)
 *  flag==1, new configuration, same number of particles
 *  flag==2, old configuration, preparing to delete (same as 0)
 *  flag==3, just inserted particle
 *  flag==4, old configuration, but only update the neighbor list
 *  flag==5, new configuration, but only update the neighbor list
 */
void Pair::updateBase(
  const vector<int> mpart,    //!< particles involved in move
  const int flag,         //!< type of move
  const char* uptype,    //!< description of update type
  vector<vector<int> > &neigh,   //!< neighbor list to update
  vector<vector<int> > &neighOne,   //!< neighbor list to update
  vector<vector<int> > &neighOneOld   //!< neighbor list to update
  ) {
  std::string uptypestr(uptype);
  if (uptypestr.compare("store") == 0) {
    if (flag == 0 || flag == 2 || flag == 4) {
      neighOne.clear();
    }
    if (flag == 0 || flag == 4) {
      neighOneOld.clear();
      if (atomCut_ == 0) {
        neighOneOld.push_back(neigh[space_->mol()[mpart[0]]]);
      } else {
        for (unsigned int impart = 0; impart < mpart.size(); ++impart) {
          for (unsigned int iii = 0; iii < neigh[mpart[impart]].size(); ++iii) {
            // cout << "adding noo " << neigh[mpart[impart]][0] << endl;
          }
          neighOneOld.push_back(neigh[mpart[impart]]);
        }
      }
    }
  }
  if (uptypestr.compare("update") == 0) {
    // update neighbor list using the precomputed neighOne vectors called during
    // energy computation
    //  this is done by replacing neigh with precomputed neighOne
    //  but also by finding differences in old and new, and updating those
    //  neighbors as well
    // cout << "updating while flag is " << flag << endl;
    if (flag == 0 || flag == 4) {
      int i = 0;
      while (i < static_cast<int>(mpart.size())) {
        int iMol = mpart[i];
        if (atomCut_ == 0) iMol = space_->mol()[iMol];
        // cout << "iMol " << iMol << endl;
        vector<int> diff;
        vector<int>::iterator it;

        // sort neighOne and neighOneOld
        std::sort(neighOne[i].begin(), neighOne[i].end());
        std::sort(neighOneOld[i].begin(), neighOneOld[i].end());

        // for all particles which are present in new neigh list that are not
        // in old neigh list
        // add those particles to jpart neighborlist
        std::set_difference(neighOne[i].begin(), neighOne[i].end(),
          neighOneOld[i].begin(), neighOneOld[i].end(), back_inserter(diff));
        for (it = diff.begin(); it != diff.end(); ++it) {
          // cout << "iMol " << iMol << " gained " << *it << endl;
          neigh[*it].push_back(iMol);
        }

        // for all particles which are present in old neigh list that are not
        // in new neigh list
        //  delete those particles from jpart neighborlist
        diff.clear();
        std::set_difference(neighOneOld[i].begin(), neighOneOld[i].end(),
          neighOne[i].begin(), neighOne[i].end(), std::back_inserter(diff));
        for (it = diff.begin(); it != diff.end(); ++it) {
          // cout << "iMol " << iMol << " lost " << *it << endl;
          for (unsigned int k = 0; k < neigh[*it].size(); ++k) {
            if (neigh[*it][k] == iMol) neigh[*it].erase(neigh[*it].begin() + k);
          }
        }

        // finally, for iMol neighbor list, replace old with new
        neigh[iMol].clear();
        for (it = neighOne[i].begin(); it != neighOne[i].end(); ++it) {
          neigh[iMol].push_back(*it);
        }

        // iterate through mpart or terminate
        if (atomCut_ == 0) {
          i = static_cast<int>(mpart.size());
        } else {
          ++i;
        }
      }
    }

    // update neighbor list using the precomputed neighOne vector called during
    // energy computation
    //  this is done by adding precomputed neighOne to neigh
    //  but also by updating neighbor lists of all other particles
    if (flag == 3) {
      int i = 0;
      while (i < static_cast<int>(mpart.size())) {
        int iMol = mpart[i];
        if (atomCut_ == 0) iMol = space_->mol()[iMol];
        vector<int>::iterator it;

        // for all neighbors of new particles, add new particles to their
        // neighbor list
        // simultaneously, add neigh list of new particles
        for (it = neighOne[i].begin(); it != neighOne[i].end(); ++it) {
          // std::cout << "iMol " << iMol << " gained " << *it << std::endl;
          neigh[*it].push_back(iMol);
          neigh[iMol].push_back(*it);
        }

        // iterate through mpart or terminate
        if (atomCut_ == 0) {
          i = static_cast<int>(mpart.size());
        } else {
          ++i;
        }
      }
    }
  }
}

/**
 *  check neigh list by re-building
 */
int Pair::checkNeigh() {
  // cout << "checking neigh of size " << neigh_.size() << endl;
  int neighMatch = 1;

  // record neigh list before rebuild in order to check that they are the same
  vector<vector<int> > neigh;
  int nMol = space_->nMol();
  if (atomCut_ == 0) {
    neigh.resize(space_->nMol(), vector<int>(0));
  } else {
    nMol = space_->natom();
    neigh.resize(space_->natom(), vector<int>(0));
  }

  // check that first dimension of neigh list is the correct size
  for (unsigned int i = 0; i < neigh_.size(); ++i) {
    if (nMol < static_cast<int>(neigh_[i].size())) {
      std::cerr << "problem " << std::endl;
      std::cerr << "this should absolutely never happen! nMol " << nMol
                << " neighsize " << neigh_[i].size() << std::endl;
      neighMatch = 0;
    }
    for (unsigned int j = 0; j < neigh_[i].size(); ++j) {
      neigh[i].push_back(neigh_[i][j]);
    }
  }

  // rebuild and check
  buildNeighList();
  if (neigh.size() != neigh_.size()) {
    std::cerr << "problem1 " << std::endl;
    std::cerr << "neighsize " << neigh.size() << " p.neighsize "
              << neigh_.size() << std::endl;
    neighMatch = 0;
  }
  for (unsigned int i = 0; i < neigh_.size(); ++i) {
    // cout << " I " << i << " size " << neigh_[i].size() << endl;
    if (neigh[i].size() != neigh_[i].size()) {
      std::cerr << "problem2 " << std::endl;
      std::cerr << "i " << i << " neighsizeOld " << neigh[i].size()
                << " neighsize " << neigh_[i].size() << std::endl;
      neighMatch = 0;
    }
    // sort on-fly neighborlist for comparison
    std::sort(neigh[i].begin(), neigh[i].end());
    for (unsigned int j = 0; j < neigh_[i].size(); ++j) {
      if (neigh[i][j] != neigh_[i][j]) {
        std::cerr << "problem3 " << std::endl;
        std::cerr << "j " << i << " neighj " << neigh[i][j] << " p.neighj "
                  << neigh_[i][j] << std::endl;
        neighMatch = 0;
      }
    }
  }
  ASSERT(neighMatch != 0, "neighbor list doesn't match");
  return neighMatch;
}

/**
 * check that energy of configuration and running energy of individual particles
 * match
 *  flag=0, check current peTot_ vs recomputed with initEnergy()
 *  flag=1, check peTot_ from initEnergy() vs peTot_ from multiPartEner with
 *  mpart=all atoms
 *  flag=2, check COM forces, fCOM_ from initEnergy() vs allPartEnerForce()
 *
 */
int Pair::checkEnergy(const double tol,   //!< tolerance for energies to match
  const int flag) {     //!< type of check
  double peTot1 = 0., peTot2 = 0.;

  // flag == 0 if checking current energy from peTot() vs recomputed energy
  // after "Forces"
  if (flag == 0) {
    peTot1 = peTot();
    initEnergy();
    peTot2 = peTot();
    ASSERT(fabs(peTot1 - peTot2) < tol,
      std::setprecision(std::numeric_limits<double>::digits10+2)
      << "checkEnergy(flag(" << flag << ")), |peTot1(" << peTot1
      << ") from running - peTot2(" << peTot2 << ") from initEnergy() = "
      << peTot1 - peTot2 << "| > tol(" << tol << ").");

  // flag == 1 if computing energy from scratch, first with Forces,
  // then multiPartEner against all molecules (half)
  } else if (flag == 1) {
    initEnergy();
    peTot1 = peTot();

    space_->xMolGen();
    for (int iMol = 0; iMol < space_->nMol(); ++iMol) {
      vector<int> mpart;
      for (int isite = space_->mol2part()[iMol];
           isite < space_->mol2part()[iMol+1]; ++isite) {
        mpart.push_back(isite);
      }
      peTot2 += 0.5*multiPartEner(mpart, 0);
    }
    if (fabs(peTot1 - peTot2) > tol) {
      ASSERT(0, std::setprecision(std::numeric_limits<double>::digits10+2)
        << "checkEnergy(flag(" << flag << ")), |peTot1(" << peTot1
        << ") from ForceS() - peTot2(" << peTot2 << ") from multiPartEner = "
        << peTot1 - peTot2 << "| > tol(" << tol << ").");
      return 0;
    }
  } else if (flag == 2) {
    // flag == 2 if computing COM forces from initEnergy() and comparing with
    // allPartEnerForce()
    std::ostringstream ermesg;
    initEnergy();
    vector<vector<double> > fForces = fCOM_;
    allPartEnerForce(1);
    bool fail = false;
    for (int iMol = 0; iMol < space_->nMol(); ++iMol) {
      for (int dim = 0; dim < space_->dimen(); ++dim) {
        if (fabs(fForces[iMol][dim] - fCOM_[iMol][dim]) > tol) {
          ermesg << "f[iMol=" << iMol << "][dim=" << dim  << "] diff is "
            << fabs(fForces[iMol][dim] - fCOM_[iMol][dim]) << " fForces "
            << fForces[iMol][dim] << " f " << fCOM_[iMol][dim] << endl;
          fail = true;
        }
      }
    }
    ASSERT(!fail, std::setprecision(std::numeric_limits<double>::digits10+2)
      << "checkEnergy(flag(" << flag << ")), forces between initEnergy() and"
      << "allPartEnerForce() do not match" << endl << ermesg.str() );
  } else {
    ASSERT(0, "unrecognized flag(" << flag << ") in checkEnergy");
  }
  return 1;
}

/**
 * Initialize with LAMMPS data file
 * Reads pair parameters
 */
void Pair::initLMPData(const string fileName) {  //!< LAMMPS Data file name
  // open LAMMPS data file
  std::ifstream file(fileName.c_str());
  ASSERT(file.good(), "cannot find lammps DATA file " << fileName);

  initDataFiles_.push_back(fileName);

  // skip all lines beginning with the character "#"
  skipCharsInFile('#', file);

  // read number of atom types
  int natype;
  std::string descript;
  file >> natype >> descript;
  while (descript.compare("atom") != 0) {
    file >> natype >> descript;
  }

  // read until Pair Coeffs
  readUntil("Pair Coeffs", file);

  // read eps, sig and (optionally) sigref
  vector<double> eps(natype), sig(natype), sigref;
  if (sigrefFlag_ == 1) sigref.resize(natype);
  string line;
  for (int i = 0; i < natype; ++i) {
    int tmp;
    file >> tmp >> eps[i] >> sig[i];
    if (sigrefFlag_ == 1) file >> sigref[i];
    getline(file, line);
  }

  // assign eps, sig
  initPairData(natype, eps, sig, sigref);
}

/**
 * initialize pair data based on number of at types, eps and sig
 */
void Pair::initPairData(const int natype,
  const vector<double> eps,
  const vector<double> sig,
  const vector<double> sigref) {
  if (natype == space_->nParticleTypes()) {
    initPairParam(eps, sig, sigref);
  } else {
    const int natypeall = space_->nParticleTypes();
    ASSERT(natype <= natypeall, "number of particle types is greater than"
      << "number of particle types in space class(" << space_->nParticleTypes()
      << ").");

    // determine size of eps, sig
    int natypeprev = 0;
    if (sigepsDefault_ == 0) natypeprev += eps_.size();
    const int ns = natypeprev+natype;
    vector<double> epsfull(ns), sigfull(ns), sigreffull(ns);

    // initialize eps,sig with prexisting values
    for (int i = 0; i < natypeprev; ++i) {
      epsfull[i] = eps_[i];
      sigfull[i] = sig_[i];
      if (sigrefFlag_ == 1) sigreffull[i] = sigRef_[i];
    }

    // add new values to eps, sig
    for (int i = natypeprev; i < natypeprev+natype; ++i) {
      // int tmp;
      // file >> tmp >> epsfull[i] >> sigfull[i];
      epsfull[i] = eps[i - natypeprev];
      sigfull[i] = sig[i - natypeprev];
      if (sigrefFlag_ == 1) sigreffull[i] = sigref[i - natypeprev];
    }

    // initialize eps, sig
    initPairParam(epsfull, sigfull, sigreffull);
  }
}

/**
 * stores, restores or updates variables to avoid recompute of entire
 * configuration after every change
 */
void Pair::update(const vector<int> mpart,    //!< particles involved in move
  const int flag,         //!< type of move
  const char* uptype    //!< description of update type
  ) {
  if (neighOn_) {
    // rebuilt neighlist if all particles are updated
    if (static_cast<int>(mpart.size()) != space_->natom()) {
      updateBase(mpart, flag, uptype, neigh_, neighOne_, neighOneOld_);
    }
  }
}

/**
 * set pair interaction parameters manually. Originally created as easy way to
 * turn off interactions
 */
void Pair::epsijset(const int i, const int j, const double eps) {
  epsij_.at(i).at(j) = eps;
  vector<double> epsijsetrec;
  epsijsetrec.push_back(i);
  epsijsetrec.push_back(j);
  epsijsetrec.push_back(eps);
  epsijsetRecord_.push_back(epsijsetrec);
}

/**
 * excluded volume of all molecules in space, using sig_
 */
double Pair::exVol(
  const int nGrid,    //!< number of grid points in each dimension
  const double dProbe                  //!< diameber of probe
  ) {
  ASSERT(dimen_ == 3, "exVol requires 3 dimensions");
  const double voldiff = space_->minl() - pow(space_->vol(),
    1./static_cast<double>(dimen_));
  ASSERT(fabs(voldiff) <= 1e-12,
    "box must be a cube to compute excluded volume");
  ASSERT(nGrid % 2 == 0, "nGrid(" << nGrid << ") must be even");

  // search through equally-spaced grid of points to find fraction that overlap
  // with molecules, based on sig
  const vector<double> &x = space_->x();
  const vector<double> l = space_->l();
  const double dGrid = l[0]/static_cast<double>(nGrid);
  const int natom = space_->natom();
  double dx, dy, dz, r2;
  int overlaps = 0;
  for (int xi = 0; xi < nGrid; ++xi) {
  for (int yi = 0; yi < nGrid; ++yi) {
  for (int zi = 0; zi < nGrid; ++zi) {
    const double xx = (xi + 0.5 - nGrid/2.) * dGrid;
    const double yy = (yi + 0.5 - nGrid/2.) * dGrid;
    const double zz = (zi + 0.5 - nGrid/2.) * dGrid;
    bool overlap = false;
    int iAtom = 0;
    while ( (!overlap) && (iAtom < natom) ) {
      dx = xx - x[dimen_*iAtom];
      dy = yy - x[dimen_*iAtom+1];
      dz = zz - x[dimen_*iAtom+2];
      r2 = dx*dx + dy*dy + dz*dz;
//      cout << "r2 " << r2 << " " ;
//      cout << "xx " << xx << " " << yy << " " << zz << " ";
//      cout << "type " << space_->type()[iAtom] << " ";
      const double sig = sig_[space_->type()[iAtom]];
      double sigtmp;
      if (fabs(sig) < doubleTolerance) {
        sigtmp = 0.;
      } else {
        sigtmp = 0.5*(sig + dProbe);
      }
      if (r2 < sigtmp*sigtmp) {
        overlap = true;
        ++overlaps;
        if (xi == 0 || yi == 0 || zi == 0 ||
            xi == nGrid -1 || yi == nGrid -1 || zi == nGrid - 1)
            ASSERT(0, "box not large enough for exVol");
      }
      ++iAtom;
    }
  }}}

  // use fraction of overlapping points to compute volume from boundary
  return space_->vol()*overlaps/pow(nGrid, dimen_);
}

/**
 * update potential energy
 *
 * simple pair updating method doesn't work with neighborlist. If you're
 * attempt AVB with ConfigBias, likely you forgot to enable dual-cut.
 * However, this may work with TrialGCA
 */
void Pair::update(const double de) {
  peTot_ += de;
}

/**
 * Write positions to a file
 */
int Pair::printxyz(const char* fileName,  //!< file with configuration
  const int initFlag,  //!< open if flag is 1, append if flag is 0
  const std::string comment) {
  // xyz file assumes 3D, but <3D is ok because you can simply define a plane
  // if floppy box, print in GRO format instead
  if (space_->floppyBox() == 1) return printGRO(fileName, initFlag);

  stringstream ss;
  ss << fileName << ".xyz";
  FILE * xyzFile = NULL;
  if (initFlag == 1) {
    fileBackUp(ss.str().c_str());
    xyzFile = fopen(ss.str().c_str(), "w");
  } else if (initFlag == 2) {
    if (fileExists(ss.str().c_str())) {
      return 0;
    } else {
      xyzFile = fopen(ss.str().c_str(), "w");
    }
  } else if (initFlag == 0) {
    xyzFile = fopen(ss.str().c_str(), "a");
  } else {
    ASSERT(0, "Unrecognized initFlag\n");
  }

  // read only from space class
  const int natom = space_->natom();
  const vector<double> x = space_->x();
  const vector<int> type = space_->type();

  bool nonInteractingSite = false;
  fprintf(xyzFile, "%i\n%f ", natom, order());
  for (int dim = 0; dim < space_->dimen(); ++dim) {
    fprintf(xyzFile, "%f ", space_->l()[dim]);
  }
  fprintf(xyzFile, "%f %s\n", space_->xyTilt(), comment.c_str());
  if (xyzFile != NULL) {
    for (int ipart = 0; ipart < natom; ++ipart) {
      if (eps_[type[ipart]] == 0) {
        fprintf(xyzFile, "H ");
        nonInteractingSite = true;
      } else if (type[ipart] == 1) {
        fprintf(xyzFile, "O ");
      } else if (type[ipart] == 2) {
        fprintf(xyzFile, "C ");
      } else if (type[ipart] == 3) {
        fprintf(xyzFile, "N ");
      } else if (type[ipart] == 4) {
        fprintf(xyzFile, "A ");
      } else if (type[ipart] == 5) {
        fprintf(xyzFile, "B ");
      } else if (type[ipart] == 6) {
        fprintf(xyzFile, "D ");
      } else {
        fprintf(xyzFile, "H ");
      }
      for (int i = 0; i < dimen_; ++i) {
        fprintf(xyzFile, "%f ", x[dimen_*ipart+i]);
        // fprintf(xyzFile, "%24.20f ", x[dimen_*ipart+i]);
      }

      // print constant plane for <3D
      for (int i = dimen_; i < 3; ++i) {
        fprintf(xyzFile, "%f ", 0.);
      }
      fprintf(xyzFile, "\n");
    }
  }
  fclose(xyzFile);

  // write vmd script to visualize in fileName appended with .vmd
  if ( (initFlag == 1) || (initFlag == 2) ) {
    ss.str("");
    ss << fileName << ".vmd";
    std::ofstream vmdf(ss.str().c_str());
    double radius;
    vmdf << "display projection Orthographic" << endl
      << "color Display Background white" << endl
      << "axes location Off" << endl;
    if (initFlag == 2) {
      radius = 10;
      vmdf << "mol load xyz " << trim("/", fileName) << ".xyz" << endl
           << "animate delete beg 0 end 0" << endl
           << "mol addfile " << trim("/", fileName) << ".xtc waitfor all"
           << endl;
    } else {
      radius = 1;
      vmdf << "topo readvarxyz " << trim("/", fileName) << ".xyz" << endl;
    }
    vmdf << "mol modstyle 0 0 VDW 1.0000000 120.000000" << endl;
    if (space_->nParticleTypes() == 1) {
      vmdf << "set sel [atomselect top \"name N\"]" << endl
           << "$sel set radius " << 0.5*radius*sig_[0] << endl
           << "$sel set mass 1" << endl;
    }
    if (space_->nParticleTypes() > 1) {
      vmdf << "set sel [atomselect top \"name O\"]" << endl;
      if (sig_.size() > 1) {
        vmdf << "$sel set radius " << 0.5*radius*sig_[1] << endl;
      } else {
        vmdf << "$sel set radius " << 0.5*radius*sig_[0] << endl;
      }
      vmdf << "$sel set mass 1" << endl;
    }
    if (space_->nParticleTypes() > 2) {
      vmdf << "set sel [atomselect top \"name C\"]" << endl
           << "$sel set radius " << 0.5*radius*sig_[2] << endl
           << "$sel set mass 1" << endl;
    }
    if (space_->nParticleTypes() > 3) {
      vmdf << "set sel [atomselect top \"name N\"]" << endl
           << "$sel set radius " << 0.5*radius*sig_[3] << endl
           << "$sel set mass 1" << endl;
      if (space_->nParticleTypes() > 4) {
        vmdf << "set sel [atomselect top \"name A\"]" << endl
             << "$sel set radius " << 0.5*radius*sig_[4] << endl
             << "$sel set mass 1" << endl;
      }
      if (space_->nParticleTypes() > 5) {
        vmdf << "set sel [atomselect top \"name B\"]" << endl
             << "$sel set radius " << 0.5*radius*sig_[5] << endl
             << "$sel set mass 1" << endl;
      }
      if (space_->nParticleTypes() > 6) {
        vmdf << "set sel [atomselect top \"name D\"]" << endl
             << "$sel set radius " << 0.5*radius*sig_[6] << endl
             << "$sel set mass 1" << endl;
      }
      vmdf << "set sel [atomselect top \"name H\"]" << endl
           << "$sel set radius " << 0.5*radius*sig_[0] << endl
           << "$sel set mass 1" << endl;
    } else {
      vmdf << "set sel [atomselect top \"name H\"]" << endl
           << "$sel set radius " << 0.5*radius*sig_[0] << endl
           << "$sel set mass 1" << endl;
    }
    if (nonInteractingSite) {
      vmdf << "set sel [atomselect top \"name H\"]" << endl
           << "$sel set radius " << 0.2*radius*sig_[0] << endl
           << "$sel set mass 1" << endl;
    }
  }
  return 0;
}

/**
 * set i-j cutoff
 */
void Pair::rCutijset(const int itype, const int jtype, const double rCut) {
  rCutij_.resize(epsij_.size(), vector<double>(epsij_.size()));
  rCutij_[itype][jtype] = rCut;
  rCutij_[jtype][itype] = rCut;
  if (static_cast<int>(rCutMax_.size()) <= itype) {
    rCutMax_.resize(itype+1);
  }
  if (static_cast<int>(rCutMax_.size()) <= jtype) {
    rCutMax_.resize(jtype+1);
  }
  rCutMax_[itype] = *std::max_element(rCutij_[itype].begin(),
                                      rCutij_[itype].end());
  rCutMax_[jtype] = *std::max_element(rCutij_[jtype].begin(),
                                      rCutij_[jtype].end());
  rCutMaxAll_ = *std::max_element(rCutMax_.begin(), rCutMax_.end());
}

/**
 * set type of atoms to keep track of neighbors
 */
void Pair::neighTypeSet(const int itype) {
  neighTypeScreen_ = 1;
  neighType_.push_back(itype);
}

/**
 * loop through all particles which interact with mpart, using atom-based
 * cut-off to compute interaction energies
 */
double Pair::multiPartEnerAtomCut(const vector<int> mpart) {
  const bool verbose = false;
  // const bool verbose = true;
  if (verbose) cout << "in Pair::multiPartEnerAtomCut(" << endl;
  // shorthand for read-only space variables
  const vector<int> type = space_->type();
  const vector<double> &x = space_->x();
  const vector<double> &l = space_->l();
  const vector<int> &mol = space_->mol();
  const double xyTilt = space_->xyTilt(), xzTilt = space_->xzTilt(),
    yzTilt = space_->yzTilt();

  // declare variables for optimization
  double r2, xi, yi, zi, dx, dy, dz, dx0, dy0;
  const double lx = l[0], ly = l[1], lz = l[2],
    halflx = lx/2., halfly = ly/2., halflz = lz/2.;
  int ipart, jpart, itype, jtype;
  double rCutij;

  // initialize neigh, neighCut and peMap
  initNeighCutPEMap(mpart);

  // loop through atoms in mpart
  for (unsigned int ii = 0; ii < mpart.size(); ++ii) {
    ipart = mpart[ii];
    const int iMol = space_->mol()[ipart];
    if (verbose) cout << "ipart " << ipart << endl;
    itype = type[ipart];
    if (eps_[itype] != 0) {   // ignore non-interacting atoms
      xi = x[dimen_*ipart];
      yi = x[dimen_*ipart+1];
      zi = x[dimen_*ipart+2];

      // obtain neighList with cellList
      if (space_->cellType() == 1) {
        if (cheapEnergy_) {
          // cheap energy switches cut-off to sigma.
          // Neighlist must be built for largest sigma
          space_->buildNeighListCellAtomCut(ipart);
        } else {
          if (rCutMax_.size() > 0) {
            if (rCutMax_[itype] <= space_->dCellMin()) {
              space_->buildNeighListCellAtomCut(ipart);
            } else {
              space_->initCellAtomCut(1);   // set neighListChosen to all atoms
            }
          } else {
            space_->initCellAtomCut(1);   // set neighListChosen to all atoms
          }
        }
      }
      const vector<int> &neigh = space_->neighListChosen();

      // loop neighboring atoms
      for (unsigned int ineigh = 0; ineigh < neigh.size(); ++ineigh) {
        jpart = neigh[ineigh];
        const int jMol = mol[jpart];
        if ( intraCheck_(ipart, jpart, iMol, jMol) &&
             ( (ipart < jpart) || (!findInList(jpart, mpart)) ) ) {
          jtype = type[jpart];
          if (eps_[jtype] != 0) {
            // separation distance with periodic boundary conditions
            dx = dx0 = xi - x[dimen_*jpart];
            dy = dy0 = yi - x[dimen_*jpart+1];
            dz = zi - x[dimen_*jpart+2];
            if (fabs(dz) > halflz) {
              if (dz < 0.) {
                dz += lz;
                dy += yzTilt;
                dx += xzTilt;
                dy0 += yzTilt;
                dx0 += xzTilt;
              } else {
                dz -= lz;
                dy -= yzTilt;
                dx -= xzTilt;
                dy0 -= yzTilt;
                dx0 -= xzTilt;
              }
            }
            if (fabs(dy0) > halfly) {
              if (dy0 < 0.) {
                dy += ly;
                dx += xyTilt;
                dx0 += xyTilt;
              } else {
                dy -= ly;
                dx -= xyTilt;
                dx0 -= xyTilt;
              }
            }
            if (fabs(dx0) > halflx) {
              if (dx0 < 0.) {
                dx += lx;
              } else {
                dx -= lx;
              }
            }
            r2 = dx*dx+dy*dy+dz*dz;

            // no interaction beyond cut-off distance
            rCutij = rCutij_[itype][jtype];
            if (verbose) cout << "ipart " << ipart << " jpart " << jpart << " r2 " << r2
              << " rCutij " << rCutij << " itype " << itype << " jtype "
              << jtype << endl;
            if (r2 < rCutij*rCutij) {
              multiPartEnerAtomCutInner(r2, itype, jtype);
              if (verbose) cout << "mpe inner peSRone " << peSRone_ << " ip "
                << ipart << " jp " << jpart << " it " << itype << " jt "
                << jtype << endl;
              if (verbose) cout << "neighOn_ " << neighOn_ << " neighAboveSq_ "
                << neighAboveSq_ << " neighBelowSq_ " << neighBelowSq_
                << " r2 " << r2 << endl;
              // store new neighbor list
              if (neighOn_) {
                if ( (r2 < neighAboveSq_) && (r2 > neighBelowSq_) ) {
                  if (neighTypeScreen_ == 1) {
                    if ( ( (neighType_[0] == itype) &&
                           (neighType_[1] == jtype) ) ||
                         ( (neighType_[1] == itype) &&
                           (neighType_[0] == jtype) ) ) {
                      neighOne_[ii].push_back(jpart);
//                      cout << "addinga neigh " << ii << " " << jpart << " "
//                        << neighOne_.size() << endl;
                    }
                  } else {
                    neighOne_[ii].push_back(jpart);
                  }
                }
              }

              // store neighCut and peMap
              storeNeighCutPEMap(jpart, ii);
            }
          }
        }
      }
    }
  }
  if (peMapOn_ == 1) peSRone_ = peSRoneAlt_;
  return peSRone_;
}

/**
 * Returns the total scalar virial of the system
 */
double Pair::vrTot() {
  double vrTotTemp = 0.;
  for (int ipart = 0; ipart < space_->natom(); ++ipart) {
    for (int i = 0; i < space_->dimen(); ++i) {
      vrTotTemp += vr_[ipart][i][i];
    }
  }
  return vrTotTemp;
}

/**
 * potential energy and forces of all particles
 */
double Pair::allPartEnerForceAtomCutNoCell() {
  // shorthand for read-only space variables
  const int natom = space_->natom();
  const vector<double> l = space_->l();
  const vector<double> x = space_->x();
  const vector<int> mol = space_->mol();
  const vector<int> type = space_->type();
  const double xyTilt = space_->xyTilt(), xzTilt = space_->xzTilt(),
    yzTilt = space_->yzTilt();

  // declare variables for optimization
  double r2, xi, yi, zi, dx, dy, dz, dx0, dy0;
  const double lx = l[0], ly = l[1], lz = l[2],
    halflx = lx/2., halfly = ly/2., halflz = lz/2.;
  int ipart, jpart, iMol, jMol, itype, jtype;

  // loop through nearest neighbor atom pairs
  for (ipart = 0; ipart < natom - 1; ++ipart) {
    itype = type[ipart];
    if (eps_[itype] != 0) {
      xi = x[dimen_*ipart];
      yi = x[dimen_*ipart+1];
      zi = x[dimen_*ipart+2];
      iMol = mol[ipart];
      for (jpart = ipart + 1; jpart < natom; ++jpart) {
        jMol = mol[jpart];
        jtype = type[jpart];
        if (intraCheck_(ipart, jpart, iMol, jMol) && (eps_[jtype] != 0))  {
          // separation distance with periodic boundary conditions
          dx = dx0 = xi - x[dimen_*jpart];
          dy = dy0 = yi - x[dimen_*jpart+1];
          dz  = zi - x[dimen_*jpart+2];
          if (fabs(dz) > halflz) {
            if (dz < 0.) {
              dz += lz;
              dy += yzTilt;
              dx += xzTilt;
              dy0 += yzTilt;
              dx0 += xzTilt;
            } else {
              dz -= lz;
              dy -= yzTilt;
              dx -= xzTilt;
              dy0 -= yzTilt;
              dx0 -= xzTilt;
            }
          }
          if (fabs(dy0) > halfly) {
            if (dy0 < 0.) {
              dy += ly;
              dx += xyTilt;
              dx0 += xyTilt;
            } else {
              dy -= ly;
              dx -= xyTilt;
              dx0 -= xyTilt;
            }
          }
          if (fabs(dx0) > halflx) {
            if (dx0 < 0.) {
              dx += lx;
            } else {
              dx -= lx;
            }
          }
          r2 = dx*dx+dy*dy+dz*dz;

          // no interaction beyond cut-off distance
          const double rCutij = rCutij_[itype][jtype];
          if (r2 < rCutij*rCutij) {
            allPartEnerForceInner(r2, dx, dy, dz, itype, jtype, iMol, jMol);
          }
        }
      }
    }
  }
  return peSRone_;
}

/**
 * potential energy and forces of all particles
 */
double Pair::allPartEnerForceAtomCutNoCell2D() {
  // shorthand for read-only space variables
  const int natom = space_->natom();
  const vector<double> l = space_->l();
  const vector<double> x = space_->x();
  const vector<int> mol = space_->mol();
  const vector<int> type = space_->type();
  const double xyTilt = space_->xyTilt();

  // declare variables for optimization
  double r2, xi, yi, dx, dy;
  const double lx = l[0], ly = l[1], halflx = lx/2., halfly = ly/2.;
  int ipart, jpart, iMol, jMol, itype, jtype;

  // loop through nearest neighbor atom pairs
  for (ipart = 0; ipart < natom - 1; ++ipart) {
    itype = type[ipart];
    if (eps_[itype] != 0) {
      xi = x[dimen_*ipart];
      yi = x[dimen_*ipart+1];
      iMol = mol[ipart];
      for (jpart = ipart + 1; jpart < natom; ++jpart) {
        jMol = mol[jpart];
        jtype = type[jpart];
        if (intraCheck_(ipart, jpart, iMol, jMol) && (eps_[jtype] != 0))  {
          // separation distance with periodic boundary conditions
          dx = xi - x[dimen_*jpart];
          dy = yi - x[dimen_*jpart+1];
          if (fabs(xyTilt) < doubleTolerance) {
            if (dx >  halflx) dx -= lx;
            if (dx < -halflx) dx += lx;
            if (dy >  halfly) dy -= ly;
            if (dy < -halfly) dy += ly;
          } else {
            if (fabs(dy) > halfly) {
              if (dy < 0.) {
                dy += ly;
                dx += xyTilt;
              } else {
                dy -= ly;
                dx -= xyTilt;
              }
            }
            if (fabs(dx) > halflx) {
              if (dx < 0.) {
                dx += lx;
              } else {
                dx -= lx;
              }
            }
          }
          r2 = dx*dx + dy*dy;

          // no interaction beyond cut-off distance
          const double rCutij = rCutij_[itype][jtype];
          if (r2 < rCutij*rCutij) {
            allPartEnerForceInner(r2, dx, dy, 0., itype, jtype, iMol, jMol);
          }
        }
      }
    }
  }
  return peSRone_;
}

/**
 * potential energy and forces of all particles
 */
double Pair::allPartEnerForceAtomCutCell() {
  // shorthand for read-only space variables
  const vector<double> l = space_->l();
  const vector<double> &x = space_->x();
  const vector<int> &mol = space_->mol();
  const vector<int> &type = space_->type();
  const vector<vector<int> > neighCell = space_->neighCell();
  const vector<vector<int> > cellList = space_->cellList();

  // declare variables for optimization
  double r2, xi, yi, zi, dx, dy, dz;
  const double lx = l[0], ly = l[1], lz = l[2],
    halflx = lx/2., halfly = ly/2., halflz = lz/2.;
  int ipart, jpart, iMol, jMol, itype, jtype;

  // loop through cells
  for (int iCell = 0; iCell < space_->nCell(); ++iCell) {
    // loop through atoms in cell
    for (unsigned int ip = 0; ip < cellList[iCell].size(); ++ip) {
      ipart = cellList[iCell][ip];
      itype = type[ipart];
      if (eps_[itype] != 0) {
        xi = x[dimen_*ipart];
        yi = x[dimen_*ipart+1];
        zi = x[dimen_*ipart+2];
        iMol = mol[ipart];

        // loop through neighboring cells with index >= current cell
        for (unsigned int j = 0; j < neighCell[iCell].size(); ++j) {
          const int jCell = neighCell[iCell][j];
          if (jCell >= iCell) {
            // loop through atoms in neighboring cell
            for (unsigned int jp = 0; jp < cellList[jCell].size(); ++jp) {
              jpart = cellList[jCell][jp];
              jMol = mol[jpart];
              jtype = type[jpart];
              if ( intraCheck_(ipart, jpart, iMol, jMol) &&
                   (eps_[jtype] != 0) &&
                   ((jCell != iCell) || (ipart < jpart) ) ) {
                // separation distance with periodic boundary conditions
                dx = xi - x[dimen_*jpart];
                dy = yi - x[dimen_*jpart+1];
                dz = zi - x[dimen_*jpart+2];
                if (dx >  halflx) dx -= lx;
                if (dx < -halflx) dx += lx;
                if (dy >  halfly) dy -= ly;
                if (dy < -halfly) dy += ly;
                if (dz >  halflz) dz -= lz;
                if (dz < -halflz) dz += lz;
                r2 = dx*dx + dy*dy + dz*dz;

                // no interaction beyond cut-off distance
                const double rCutij = rCutij_[itype][jtype];
                if (r2 < rCutij*rCutij) {
                  allPartEnerForceInner(r2, dx, dy, dz, itype, jtype,
                                        iMol, jMol);
                }
              }
            }
          }
        }
      }
    }
  }
  return peSRone_;
}

/**
 * loop through all particles which interact with mpart, using atom-based
 * cut-off to compute interaction energies
 */
double Pair::multiPartEnerAtomCut2D(const vector<int> mpart) {
  // cout << "in multiPartEnerAtomCut2D" << endl;
  // shorthand for read-only space variables
  const vector<double> &x = space_->x();
  const vector<double> &l = space_->l();
  const vector<int> &mol = space_->mol();
  const vector<int> &type = space_->type();
  const double xyTilt = space_->xyTilt();

  // declare variables for optimization
  double r2, xi, yi, dx, dy;
  const double lx = l[0], ly = l[1], halflx = lx/2., halfly = ly/2.;

  // initialize neigh, neighCut and peMap
  initNeighCutPEMap(mpart);

  // loop through all particles in mpart
  for (unsigned int ii = 0; ii < mpart.size(); ++ii) {
    const int ipart = mpart[ii];
    if (eps_[type[ipart]] != 0) {
      const int iMol = mol[ipart];
      const int iType = type[ipart];
      xi = x[dimen_*ipart];
      yi = x[dimen_*ipart+1];

      // obtain neighList with cellList
      if (space_->cellType() == 1) {
        if (cheapEnergy_) {
          // cheap energy switches cut-off to sigma. Neighlist must be built
          // for largest sigma
          space_->buildNeighListCellAtomCut(ipart);
        } else {
          if (rCutMax_.size() > 0) {
            if (rCutMax_[iType] <= space_->dCellMin()) {
              space_->buildNeighListCellAtomCut(ipart);
            } else {
              space_->initCellAtomCut(1);   // set neighListChosen to all atoms
            }
          } else {
            space_->initCellAtomCut(1);   // set neighListChosen to all atoms
          }
        }
      }
      const vector<int> &neigh = space_->neighListChosen();

      // loop neighboring atoms
      for (unsigned int ineigh = 0; ineigh < neigh.size(); ++ineigh) {
        const int jpart = neigh[ineigh];
        const int jType = type[jpart];
        if ( intraCheck_(ipart, jpart, iMol, mol[jpart]) &&
             (eps_[jType] != 0) ) {
          // separation distance with periodic boundary conditions
          dx = xi - x[dimen_*jpart];
          dy = yi - x[dimen_*jpart+1];
          if (fabs(xyTilt) < doubleTolerance) {
            if (dx >  halflx) dx -= lx;
            if (dx < -halflx) dx += lx;
            if (dy >  halfly) dy -= ly;
            if (dy < -halfly) dy += ly;
          } else {
            if (fabs(dy) > halfly) {
              if (dy < 0.) {
                dy += ly;
                dx += xyTilt;
              } else {
                dy -= ly;
                dx -= xyTilt;
              }
            }
            if (fabs(dx) > halflx) {
              if (dx < 0.) {
                dx += lx;
              } else {
                dx -= lx;
              }
            }
          }
          r2 = dx*dx + dy*dy;
          // cout << "ipart " << ipart << " jpart " << jpart << " r2 " << r2
          //  << " rcutsq " << rCutSq_ << endl;

          // no interaction beyond cut-off distance
          const double rCutij = rCutij_[iType][jType];
          if (r2 < rCutij*rCutij) {
            multiPartEnerAtomCutInner(r2, iType, jType);
          }

          // store neighCut and peMap
          storeNeighCutPEMap(jpart, ii);
        }
      }
    }
  }
  return peSRone_;
}

/**
 * return list of molecules in neighCutOne and their total potential energy
 * interactions. Note that peMap contains the cumulative peSRone_
 */
void Pair::neighCutMolPEMap(vector<int> &neigh, vector<double> &peMap) {
  ASSERT( (atomCut_ == 1) || (className_ == "PairPatchKF"),
    "neighCutMolPEMap only implemented for atomCut(" << atomCut_ << ")==1");
  ASSERT(neighCutOne_.size() == peMap_.size(),
    "non matching vec sizes, neighCutOne(" << neighCutOne_.size()
    << "), peMap(" << peMap_.size() << " in neighCutMolPEMap");

//  double peOld = 0;
  for (unsigned int ii = 0; ii < neighCutOne_.size(); ++ii) {
    for (unsigned int jj = 0; jj < neighCutOne_[ii].size(); ++jj) {
      const int iMol = space_->mol()[neighCutOne_[ii][jj]];

      // accumulated peSRone_ is stored in peMap, so subtract old value
      // const double pe = peMap_[ii][jj] - peOld;
      const double pe = peMap_[ii][jj];
      // cout << "pemap " << peMap_[ii][jj] << " iMol " << iMol << endl;

      // if iMol already listed in neigh, add potential energy
      int index = -1;
      if (findInList(iMol, neigh, index)) {
        peMap[index] += pe;
        // cout << "add pe " << pe << " of iMol" << iMol << endl;

      // else, pushback
      } else {
        // cout << "new iMol" << iMol << " w/pe " << pe << endl;
        neigh.push_back(iMol);
        peMap.push_back(pe);
      }
//      peOld = peMap_[ii][jj];
    }
  }
  ASSERT(neigh.size() == peMap.size(), "non matching mol vec sizes, neigh("
    << neigh.size() << "), peMap(" << peMap.size() << " in neighCutMolPEMap");
}

/*
 *  initialize neighCut and peMap
 */
void Pair::initNeighCutPEMap(const vector<int> mpart) {
  if (neighOn_) {
    neighOne_.clear();
    neighOne_.resize(static_cast<int>(mpart.size()), vector<int>());
  }
  if (neighCutOn_ == 1) {
    neighCutOne_.clear();
    neighCutOne_.resize(static_cast<int>(mpart.size()), vector<int>());
  }
  if (peMapOn_ == 1) {
    peSRoneAlt_ = 0;
    peMap_.clear();
    peMap_.resize(static_cast<int>(mpart.size()), vector<double>());
  }
}

/**
 * store neighCut and peMap
 */
void Pair::storeNeighCutPEMap(const int jpart,    //!< neighbor index
                              const int ii) {     //!< index of site in molecule
  if (neighCutOn_ == 1) {
    neighCutOne_[ii].push_back(jpart);
  }
  if (peMapOn_ == 1) {
    // cout << "storing peSRone_ " << peSRone_ << endl;
    peMap_[ii].push_back(peSRone_);
    peSRoneAlt_ += peSRone_;
    peSRone_ = 0;
  }
}

/**
 * set the order parameter
 */
void Pair::setOrder(const double order) {
  // initialize order parameter, subject to bounds
  orderOn_ = 1;
  order_ = order;
  if (order_ > orderMax_) {
    order_ = orderMax_;
  } else if (order_ < orderMin_) {
    order_ = orderMin_;
  }

  // check for special order parameter names
  if (orderName_.compare("MMsigma3") == 0) {
    sig_[2] = order;
    initPairParam(eps_, sig_, sigRef_);
  } else if (orderName_.compare("angle0") == 0) {
    space_->modBondAngle(0, order, space_->addMolListType().front().c_str());
//    space_->setAtomTypeAsCOM(0);
  } else if (orderName_.compare("MMsigSheetVes") == 0) {
    sig_[1] = 0.9 + 0.2*order;
    sig_[2] = 0.9 - 0.2*order;
    initPairParam(eps_, sig_, sigRef_);
  }
}

/**
 * read name of molecule as order parameter
 *  http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/xyzplugin.html
 *  [ # optional comment line ] comment line (can be blank)
    [ N                       ] # of atoms, required by this xyz reader
    [ molecule name           ] name of molecule (can be blank)
    atom1 x y z [optional data] atom name followed by xyz coords
    atom2 x y z [ ...         ] and and (optionally) other data.
 */
void Pair::readXYZOrder(std::ifstream& file) {  //!< XYZ file
  // read first two lines
  int iAtom;
  string line;
  getline(file, line);
  if (!file.eof()) {
    file >> iAtom;
    getline(file, line);
    setOrder(stod(line));

    // skip coordinates
    for (int i = 0; i < iAtom; ++i) getline(file, line);
  }
}

/**
 * compute the average number of neighbors within cutoff using peMap
 */
double Pair::avNumNeighCut() {
  // for each molecule, compute the pemap to get num neighs
  Accumulator neighs;
  for (int iMol = 0; iMol < space_->nMol(); ++iMol) {
    const vector<int> neigh = iMol2neigh(iMol);
    neighs.accumulate(static_cast<int>(neigh.size()));
  }
  return neighs.average();
}

/**
 * return a vector list of neighboring molecules of iMol, based on pair rCut
 */
vector<int> Pair::iMol2neigh(const int iMol) {
  initNeighCut(1);
  initPEMap(1);
  vector<int> mpart = space_->imol2mpart(iMol);
  multiPartEner(mpart, 0);
  vector<int> neigh;     // list of interacting neighbors
  vector<double> peMap;  // list of energy interactions in same order as neigh
  neighCutMolPEMap(neigh, peMap);
  initNeighCut(0);
  initPEMap(0);
  return neigh;
}

/**
 * return a list of minimum angles between neighbors with iMol as vertex, based on pair rCut
 */
vector<double> Pair::iMol2neighAngles(const int iMol) {
  ASSERT(space_->dimen() == 2,
    "iMol2neighAngles hardcoded for dim(" << space_->dimen() << ")=2");
  const vector<int> neigh = iMol2neigh(iMol);

  vector<double> angles;
  for (int i = 0; i < static_cast<int>(neigh.size())-1; ++i) {
    const int ineigh = neigh[i];
    for (unsigned int j = i+1; j < neigh.size(); ++j) {
      const int jneigh = neigh[j];

      // cout << "i " << i << " j " << j << endl;
      // find the vector from iMol to each of the neighbors
      vector<double> v1(space_->dimen()), v2 = v1;
      for (int dim = 0; dim < space_->dimen(); ++dim) {
        v1[dim] = space_->x(space_->mol2part()[ineigh], dim)
                - space_->x(space_->mol2part()[iMol], dim);
        v2[dim] = space_->x(space_->mol2part()[jneigh], dim)
                - space_->x(space_->mol2part()[iMol], dim);
      }
      vector<double> dx1 = space_->pbc(v1),
                     dx2 = space_->pbc(v2);
      for (int dim = 0; dim < space_->dimen(); ++dim) {
        v1[dim] += dx1[dim];
        v2[dim] += dx2[dim];
      }
      const double v1l = sqrt(vecDotProd(v1, v1)),
                   v2l = sqrt(vecDotProd(v2, v2));

      // obtain the angle
      const double angle = acos(vecDotProd(v1, v2)/v1l/v2l);
      // cout << "angle(deg) " << angle/PI*180 << endl;
      angles.push_back(angle);
    }
  }

  // sort and retain only the smallest angles
  std::sort(angles.begin(), angles.end());
  // cout << "sizes " << angles.size() << " " << neigh.size() << endl;
  const int npops = static_cast<int>(angles.size() - neigh.size());
  for (int i = 0; i < npops; ++i) angles.pop_back();
  return angles;
}

/*
 * compute intramolecular interactions
 * exclude interactions between pairs of sites that share a bond length
 *  parameter
 */
double Pair::peIntra(const int iAtom) {
  ASSERT(0, "Pair::peIntra is depreciated");
  const vector<int> &mol2part = space_->mol2part();
  const int iType = space_->type()[iAtom];
  const int iMol = space_->mol()[iAtom];
  const vector<double> &x = space_->x();
  const vector<double> &l = space_->l();

  const int firstAtom = mol2part[iMol];

  // list bonds involving iAtom
  vector<vector<int> > bList = space_->listBonds(iAtom);
  cout << "here we go: " << peSRone_ << endl;

  for (int jAtom = mol2part[iMol]; jAtom < mol2part[iMol+1]; ++jAtom) {
    if (jAtom != iAtom) {
      // skip if there is a bond between iAtom and jAtom
      int bond = 0;
      for (unsigned int iBond = 0; iBond < bList.size(); ++iBond) {
        const int b1Atom = bList[iBond][1] + firstAtom;
        const int b2Atom = bList[iBond][2] + firstAtom;
        cout << "i " << iAtom << " j " << jAtom << " b1 " << b1Atom << " b2 "
             << b2Atom << endl;
        if ( ( (b1Atom == iAtom) && (b2Atom == jAtom) ) ||
             ( (b1Atom == jAtom) && (b2Atom == iAtom) ) ) {
          bond = 1;
        }
      }
      if (bond == 0) {
        const int jType = space_->type()[jAtom];

        // separation vector, xij with periodic boundary conditions
        double r, r2 = 0;
        vector<double> xij(dimen_);
        for (int dim = 0; dim < dimen_; ++dim) {
          r = x[dimen_*iAtom+dim] - x[dimen_*jAtom+dim];
          if (r >  0.5 * l[dim]) r -= l[dim];
          if (r < -0.5 * l[dim]) r += l[dim];
          xij[dim] = r;
          r2 += r*r;
        }
        multiPartEnerAtomCutInner(r2, iType, jType);
        cout << "peSRone_ " << peSRone_ << " r2 " << r2 << " iAtom " << iAtom
             << " jAtom " << jAtom << endl;
      }
    }
  }
  const double pe = peSRone_;
  peSRone_ = 0;
  cout << "all done pe " << pe << endl;
  return pe;
}

/**
 * print positions in GRO file format
 *  open for first time if flag is 1
 *  append if flag is 0
 */
int Pair::printGRO(const char* fileName,  //!< name of file to print
  const int append) {  //!< flag for append/open, header and vmd scirpt options
  // open file (optional - append)
  std::stringstream ss;
  ss << fileName << ".gro";
//  assert ( (append >=0) && (append <= 3) );
  FILE * file = NULL;
  if ( (append == 0) ) {
    file = fopen(ss.str().c_str(), "a");
  } else if (append == 1) {
    file = fopen(ss.str().c_str(), "w");
  }

  // print header (optional)
  fprintf(file, "GRO file, t= 0.0\n%5i\n", space_->natom());

  // print vmd script (optional)
  if ( (append == 1) ) {
    ss.str("");
    ss << fileName << ".vmd";
    std::ofstream vmdf(ss.str());
    vmdf << "display projection Orthographic" << endl
         << "color Display Background white" << endl
         << "axes location Off" << endl
         << "mol new " << trim("/", fileName)
         << ".gro type {gro} first 0 last -1 step 1 waitfor all" << endl
         << "mol modstyle 0 0 VDW 1.0000000 10.000000" << endl;

    for (int itype = 0; itype < space_->nParticleTypes(); ++itype) {
    vmdf << "set sel [atomselect top \"name " << itype << "\"]" << endl
         << "$sel set radius " << sig_[itype]/2. << endl;
    }
    vmdf << "pbc box -center origin" << endl;
  }

  for (int iAtom = 0; iAtom < space_->natom(); ++iAtom) {
    const int iMol = space_->mol()[iAtom];
    string mt(space_->moltype()[iMol]);
    stringstream it;
    it << space_->type()[iAtom];
    if (space_->dimen() == 3) {
      fprintf(file, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", iMol, it.str().c_str(),
        it.str().c_str(), iAtom, space_->x(iAtom, 0)/10, space_->x(iAtom, 1)/10,
        space_->x(iAtom, 2)/10);
    } else {
      fprintf(file, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", iMol, it.str().c_str(),
        it.str().c_str(), iAtom, space_->x(iAtom, 0)/10, space_->x(iAtom, 1)/10,
        0.);
    }
  }
  if (space_->dimen() == 3) {
    fprintf(file, "%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n",
      space_->l(0)/10, space_->l(1)/10, space_->l(2)/10, 0., 0.,
      space_->xyTilt()/10, 0., space_->xzTilt()/10, space_->yzTilt()/10);
  } else {
    fprintf(file, "%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n",
      space_->l(0)/10, space_->l(1)/10, 0.1, 0., 0., space_->xyTilt()/10, 0.,
      0., 0.);
  }
  fclose(file);
  return 0;
}

/**
 * print a table file based on the potential
 */
int Pair::printTable(const char* tableFile,     //!< base file name of table
  const int nElements,  //!< number of table elements
  const double sigFac) {  //!< sigFac*sig = minimum separation distance in table
  // check if a tabular style is already in use
  if ( (className_.compare("PairTabular") == 0) ||
       (className_.compare("PairTabular1D") == 0)  ) {
    return 0;
  }

  // loop through each unique pair of particle types
  for (int iType = 0; iType < space_->nParticleTypes(); ++iType) {
    for (int jType = iType; jType < space_->nParticleTypes(); ++jType) {
      // initialize table file
      stringstream ss;
      ss << tableFile << "i" << iType << "j" << jType;
      std::ofstream file(ss.str().c_str());

      // determine table spacing
      //  in this case, tabulate in terms of r^2
      const double rMin = sigFac*sigij_[iType][jType],
        rMax = rCutij_[iType][jType],
        sMin = rMin*rMin,
        sMax = rMax*rMax,
        ds = (sMax-sMin)/static_cast<double>(nElements);

      // output header
      file << std::setprecision(std::numeric_limits<double>::digits10+2)
           << "# tableType " << className_ << endl
           << "# tabDims 1" << endl
           << "# dimorder0 0" << endl
           << "# dimn0 " << nElements + 1 << endl
           << "# dimmin0 " << sMin << endl
           << "# dimmax0 " << sMax << endl;

      ASSERT(space_->nMol() == 0,
        "assume that no molecules are present to make table");

      // add first molecule on origin
      vector<double> xAdd(space_->dimen());
      space_->xAdd = xAdd;
      space_->addMol(iType);

      // second second molecule on origin
      space_->xAdd = xAdd;
      space_->addMol(jType);
      addPart();
      const vector<int> jpart = space_->imol2mpart(1);
      space_->xStore(jpart);

      // loop through separation distances
      for (int step = 0; step < nElements + 1; ++step) {
        const double s = sMin + static_cast<double>(step)*ds,
          r = sqrt(s);
        xAdd[0] = r;
        space_->transMol(1, xAdd);
        initEnergy();
        file << std::setprecision(std::numeric_limits<double>::digits10+2)
             << peTot() << endl;
        // file << peTot() << " r " << r << " s " << s << endl;
        space_->restore(jpart);
      }

      // remove the molecules
      for (int imol = 0; imol < 2; ++imol) {
        vector<int> ipart = space_->imol2mpart(0);
        delPart(ipart);
        space_->delPart(ipart);
      }
    }
  }

  return 1;
}

/**
 * set all rCutij to sigmaij
 */
void Pair::sig2rCut() {
  for (unsigned int iType = 0; iType < sigij_.size(); ++iType) {
    for (unsigned int jType = iType; jType < sigij_.size(); ++jType) {
      rCutijset(iType, jType, sigij_[iType][jType]);
    }
  }
}

#ifdef JSON_
/**
 * Initialize with JSON data file
 * Reads pair parameters
 */
void Pair::initJSONData(const string fileName) {  //!< LAMMPS Data file name
  // open a JSON file
  std::ifstream file(fileName.c_str());
  ASSERT(file.good(), "cannot find json DATA file " << fileName);
  nlohmann::json j;
  file >> j;

  initDataFiles_.push_back(fileName);

  // determine number of atom types
  Space stmp(dimen_, 0);
  stmp.initData(fileName);
  const int natype = stmp.nParticleTypes();

  // read eps, sig
  vector<double> eps(natype), sig(natype), sigref;
  if (sigrefFlag_ == 1) sigref.resize(natype);
  for (int i = 0; i < natype; ++i) {
    eps[i] = j["pairCoeffs"][i]["eps"];
    sig[i] = j["pairCoeffs"][i]["sig"];
    if (sigrefFlag_ == 1) sigref[i] = j["pairCoeffs"][i]["sigref"];
  }

  // assign eps, sig
  initPairData(natype, eps, sig, sigref);
}
#endif  // JSON_

/**
 * Initialize with data file
 */
void Pair::initData(const std::string fileName) {  //!< Data file name
  // use file extension to determine whether to use JSON or LMP data files
  if (trim(".", fileName) == "json") {
    #ifdef JSON_
      initJSONData(fileName);
    #else
      ASSERT(0, "JSON files are not implemented, but file " << fileName <<
        " was attempted to be read");
    #endif  // JSON_

  } else {
    initLMPData(fileName);
  }
}

/**
 * check if there is an intramolecular interaction that is allowed
 */
bool Pair::intraCheck_(const int ipart, const int jpart,
                       const int iMol, const int jMol) {
  if (nonphys_[jpart] != 0) return false;
  if (intra_ != 2) {
    if (jMol != iMol) return true;    // intermolecular interactions
  } else {
    if (jMol != iMol) return false;   // no inter for intra_ == 2
  }
//  cout << "checking, intra " << intra_ << endl;
  if (intra_ > 0) {
    //if (ipart > jpart) return false;
    const int iAtomFirst = space_->mol2part()[iMol];
    const int ia = ipart - iAtomFirst;
    const int ja = jpart - iAtomFirst;
    //if (space_->intraMap()[ia][ja] == 1) {
    const int molid = space_->molid()[iMol];
//    cout << "molid " << molid << endl;
//    cout << "iaf " << iAtomFirst << " ia " << ia << " ja " << ja << "map "
//         << space_->intraMap()[molid][ia][ja] << endl;
    if (space_->intraMap()[molid][ia][ja] == 1) {
//      cout << "here i am" << endl;
      return true;
    }
  }
  return false;
}

/**
 * initialize intramolecular interactions
 */
void Pair::initIntra(const int flag, vector<vector<int> > map) {
  initIntra(flag);
  space_->initIntra(map);
}

void Pair::initIntraBonded(const int flag) {
  shared_ptr<Space> space = space_->addMolList()[0];
  vector<vector<int> > bondList = space->bondList();
  vector<vector<int> > intraMap(space->natom(), vector<int>(space->natom(), 1));
  for (int iatom = 0; iatom < space->natom(); ++iatom) {
    intraMap[iatom][iatom] = 0;
  }
  for (int ibond = 0; ibond < static_cast<int>(bondList.size()); ++ibond) {
    const int iAtom = bondList[ibond][1];
    const int jAtom = bondList[ibond][2];
    intraMap[iAtom][jAtom] = 0;
    intraMap[jAtom][iAtom] = 0;
  }
  initIntra(flag, intraMap);
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_



