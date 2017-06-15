#include "./pair_tabular_1d.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

PairTabular1D::PairTabular1D(Space* space)
  : Pair(space, 0.) {
  defaultConstruction();
}
PairTabular1D::PairTabular1D(Space* space,
  const char* fileName)
  : Pair(space, fileName) {
  defaultConstruction();
  string tabfile = fstos("tabFileName", fileName);
  readTable(tabfile.c_str());
}

/**
 * defaults in constructor
 */
void PairTabular1D::defaultConstruction() {
  className_.assign("PairTabular1D");
  tol_ = 0.75;
}

/**
 * write restart file
 */
void PairTabular1D::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  file << "# tabFileName " << tabFileName_ << endl;
}

/**
 */
int PairTabular1D::initEnergy() {
  // shorthand for read-only space variables
  const int natom = space_->natom();
  const vector<double> l = space_->l();
  const vector<double> x = space_->x();
  const vector<int> mol = space_->mol();
  const vector<int> type = space_->type();

  // zero accumulators: potential energy, force, and virial
  deSR_ = 0;
  fill(0., f_);
  //fCOM_.clear();
  //fCOM_.resize(space_->nMol(), vector<double>(dimen_, 0.));

  // loop through nearest neighbor atom pairs
  for (int ipart = 0; ipart < natom - 1; ++ipart) {
    if (eps_[type[ipart]] != 0) {
      for (int jpart = ipart + 1; jpart < natom; ++jpart) {
        if ( ( (intra_ != 0) || (mol[ipart] != mol[jpart]) ) &&
             (eps_[type[jpart]] != 0) ) {
          // separation vector, xij with periodic boundary conditions
          vector<double> xij(dimen_);
          for (int i = 0; i < dimen_; ++i) {
            xij[i] = x[dimen_*ipart+i] - x[dimen_*jpart+i];
          }
          const vector<double> dx = space_->pbc(xij);
          for (int dim = 0; dim < dimen_; ++dim) {
            xij[dim] += dx[dim];
          }
          double r2 = vecDotProd(xij, xij);

          // cout << "r2 " << r2 << " i " << ipart << " j " << jpart
          //   << " itype " << type[ipart] << " jtype " << type[jpart] << endl;

          // no interaction beyond cut-off distance
          if (r2 < pow(rCutij_[type[ipart]][type[jpart]], 2.)) {
            if (r2 < pow(rCutInner_[type[ipart]][type[jpart]], 2.)) {
              deSR_ += NUM_INF;
              ASSERT(0, "inner table value reached");
            } else {
              // obtain force from the potential energy table if using spline
              if (peTable_[type[ipart]][type[jpart]]->
                interpolator().compare("gslspline") == 0) {
                double fij = 0.;
                const double pe =
                  peTable_[type[ipart]][type[jpart]]->interpolate(r2, &fij);
                deSR_ += pe;
                // fij returns the derivative with respect to r2
                fij *= -2;
                for (int dim = 0; dim < dimen_; ++dim) {
                  f_[ipart][dim] += fij*xij[dim];
                  f_[jpart][dim] -= fij*xij[dim];
                }
              } else {
                const double pe =
                  peTable_[type[ipart]][type[jpart]]->interpolate(r2);
                deSR_ += pe;
              }
            }
          }
        }
      }
    }
  }
  peTot_ = deSR_;
  return 0;
}

/**
 * potential energy contribution due to particles
 */
double PairTabular1D::multiPartEner(
  const vector<int> mpart,    //!< particles to calculate energy interactions
  const int flag) {     //!< place holder for other pair styles
  if (flag == 0) {}  // remove unused parameter warning

  // zero potential energy contribution of particle ipart
  peSRone_ = 0;

  if (dimen_ == 2) {
    return multiPartEnerAtomCut2D(mpart);
  } else if (dimen_ == 3) {
    return multiPartEnerAtomCut(mpart);
  } else {
    ASSERT(0, "dimen(" << dimen_ << ") not recognized");
  }
  return 0;
}

/**
 * inner loop
 */
void PairTabular1D::multiPartEnerAtomCutInner(
  const double &r2, const int &itype, const int &jtype) {
  if (r2 < pow(rCutInner_[itype][jtype], 2.)) {
    peSRone_ += NUM_INF;
  } else {
    peSRone_ += peTable_[itype][jtype]->interpolate(r2);
  }
}

/**
 * inner loop for potential energy and forces of all particles
 */
void PairTabular1D::allPartEnerForceInner(const double &r2, const double &dx,
  const double &dy, const double &dz, const int &itype, const int &jtype,
  const int &iMol, const int &jMol) {
  if (dx*dy*dz*iMol*jMol == 0) {}  // remove unused variable warning
  multiPartEnerAtomCutInner(r2, itype, jtype);
}

/**
 * stores, restores or updates variables to avoid recompute of entire
 * configuration after every change
 */
void PairTabular1D::update(
  const vector<int> mpart,    //!< particles involved in move
  const int flag,         //!< type of move
  const char* uptype) {   //!< description of update type
  if (neighOn_) {
    updateBase(mpart, flag, uptype, neigh_, neighOne_, neighOneOld_);
  }
  std::string uptypestr(uptype);

  if (uptypestr.compare("store") == 0) {
    if (flag == 0 || flag == 2 || flag == 3) {
      deSR_ = peSRone_;
    }
  }

  if (uptypestr.compare("update") == 0) {
    if (flag == 0) {
      peTot_ += peSRone_ - deSR_;
    }
    if (flag == 2) {
      peTot_ -= deSR_;
    }
    if (flag == 3) {
      peTot_ += deSR_;
    }
  }
}

/**
 * potential energy and forces of all particles
 *  if flag == 0, dummy calculation
 *  if flag == 1, all config calculation
 */
double PairTabular1D::allPartEnerForce(const int flag) {
  peSRone_ = 0.;
  if (flag == 0) {
    peSRone_ = peTot();
    return peSRone_;
  } else if (dimen_ == 2) {
    return allPartEnerForceAtomCutNoCell2D();
  } else if (dimen_ == 3) {
    return allPartEnerForceAtomCutNoCell();
  }
  return 1e300;
}

/**
 * read tables
 */
void PairTabular1D::readTable(const char* fileName) {
  tabFileName_.assign(fileName);

  // resize rCutInner and tabij
  rCutInner_.resize(space_->nParticleTypes(), vector<double>(
    space_->nParticleTypes()));
  peTable_.resize(space_->nParticleTypes(), vector<shared_ptr<Table> >(
    space_->nParticleTypes()));

  // loop through each unique pair of particle types
  for (int iType = 0; iType < space_->nParticleTypes(); ++iType) {
    for (int jType = iType; jType < space_->nParticleTypes(); ++jType) {
      // initialize table file name
      stringstream ss;
      ss << fileName << "i" << iType << "j" << jType;

      // read rCutInner and rCut
      rCutInner_[iType][jType] = sqrt(fstod("dimmin0", ss.str().c_str()));
      rCutInner_[jType][iType] = rCutInner_[iType][jType];
      rCutijset(iType, jType, sqrt(fstod("dimmax0", ss.str().c_str())));

      // make the tables
      peTable_[iType][jType] = make_shared<Table>(ss.str().c_str());
      peTable_[jType][iType] = peTable_[iType][jType];
    }
  }
}

/**
 * set the interpolator
 */
void PairTabular1D::setInterpolator(const char* name) {
  for (int itype = 0; itype < static_cast<int>(peTable_.size()); ++itype) {
    for (int jtype = 0; jtype < static_cast<int>(peTable_[itype].size());
         ++jtype) {
      peTable_[itype][jtype]->setInterpolator(name);
    }
  }
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

