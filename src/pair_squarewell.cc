#include "./pair_squarewell.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

PairSquareWell::PairSquareWell(Space* space,
  const double rCut)
  : Pair(space, rCut) {
  defaultConstruction_();
}

PairSquareWell::PairSquareWell(Space* space,
  const char* fileName
  )
    : Pair(space, fileName) {
  defaultConstruction_();
}

void PairSquareWell::defaultConstruction_() {
  className_.assign("PairSquareWell");
}

void PairSquareWell::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
}

void PairSquareWell::initEnergy() {
  // cout << "in forces" << endl;

  // shorthand for read-only space variables
  const int natom = space_->natom();
  const vector<double> l = space_->l();
  const vector<double> x = space_->x();
  const vector<int> mol = space_->mol();
  const vector<int> type = space_->type();

  // zero accumulators: potential energy, force, and virial
  std::fill(pe_.begin(), pe_.end(), 0.);
  fill(0., f_);
  fill(0., vr_);
  peTot_ = 0;

  // loop through nearest neighbor atom pairs
  for (int ipart = 0; ipart < natom - 1; ++ipart) {
    if (eps_[type[ipart]] != 0) {
      for (int jpart = ipart + 1; jpart < natom; ++jpart) {
        if ( (mol[ipart] != mol[jpart]) && (eps_[type[jpart]] != 0) ) {
          // separation distance with periodic boundary conditions
          double r, r2 = 0;
          vector<double> xij(dimen_);
          xij[0] = x[dimen_*ipart  ] - x[dimen_*jpart];
          xij[1] = x[dimen_*ipart+1] - x[dimen_*jpart+1];
          xij[2] = x[dimen_*ipart+2] - x[dimen_*jpart+2];
          vector<double> dxpbc = space_->pbc(xij);
          for (int dim = 0; dim < dimen_; ++dim) {
            xij[dim] += dxpbc[dim];
          }
          r2 = vecDotProd(xij, xij);
          r = sqrt(r2);

          // no interaction beyond cut-off distance
          if (r < rCutij_[type[ipart]][type[jpart]]) {
            if (r < sigij_[type[ipart]][type[jpart]]) {
              // hard sphere if less than sigma
              peTot_ += NUM_INF;
            } else {
              // otherwise, square well
              peTot_ -= epsij_[type[ipart]][type[jpart]];
            }
          }
        }
      }
    }
  }
}

double PairSquareWell::multiPartEner(
  const vector<int> mpart,
  const int flag) {
  if (flag == 0) {}  // remove unused parameter warning

  ASSERT(atomCut_ == 1, "assumes atomCut(" << atomCut_ << ") is 1");

  // zero potential energy contribution of particle ipart
  peSRone_ = 0;

  // shorthand for read-only space variables
  const vector<int> type = space_->type();
  const vector<double> &x = space_->x();
  const vector<int> &mol = space_->mol();

  // declare variables for optimization
  double xi, yi, zi;

  // update neighbor list
  if (neighOn_) {
    neighOne_.clear();
    neighOne_.resize(static_cast<int>(mpart.size()), vector<int>());
  }

  // loop through all particles in mpart
  for (unsigned int ii = 0; ii < mpart.size(); ++ii) {
    const int ipart = mpart[ii];
    const int itype = type[ipart];
    if (eps_[itype] != 0) {   // ignore non-interacting atoms
      const int iMol = mol[ipart];
      xi = x[dimen_*ipart];
      yi = x[dimen_*ipart+1];
      zi = x[dimen_*ipart+2];

      /// obtain neighList with cellList
      if (space_->cellType() == 1) {
        space_->buildNeighListCellAtomCut(ipart);
      }
      const vector<int> &neigh = space_->neighListChosen();

      // loop neighboring atoms
      for (unsigned int ineigh = 0; ineigh < neigh.size(); ++ineigh) {
        const int jpart = neigh[ineigh];
        if (mol[jpart] != iMol) {
          const int jtype = type[jpart];
          if (eps_[jtype] != 0) {
            // separation distance with periodic boundary conditions
            vector<double> xij(dimen_);
            xij[0] = xi - x[dimen_*jpart];
            xij[1] = yi - x[dimen_*jpart+1];
            xij[2] = zi - x[dimen_*jpart+2];
            vector<double> dxpbc = space_->pbc(xij);
            for (int dim = 0; dim < dimen_; ++dim) {
              xij[dim] += dxpbc[dim];
            }
            double r2 = vecDotProd(xij, xij);

            // no interaction beyond cut-off distance
            const double rCutij = rCutij_[itype][jtype];
            if (r2 < rCutij*rCutij) {
              const double sigij = sigij_[itype][jtype];
              if (r2 < sigij*sigij) {
                // hard sphere if less than sigma
                peSRone_ += NUM_INF;
                // cout << "hs " << peSRone_ << endl;
              } else {
                // otherwise, square well
                peSRone_ -= epsij_[itype][jtype];

                // store new neighbor list
                if (neighOn_) {
                  if ( (r2 < neighAboveSq_) && (r2 > neighBelowSq_) ) {
                    if (neighTypeScreen_ == 1) {
                      if ( (neighType_.front() == itype) && (
                           neighType_.front() == jtype) ) {
                        neighOne_[ii].push_back(jpart);
                      }
                    } else {
                      neighOne_[ii].push_back(jpart);
                    }
                  }
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

void PairSquareWell::update(
  const vector<int> mpart,
  const int flag,
  const char* uptype
  ) {
  if (neighOn_) {
    // rebuilt neighlist if all particles are updated
    if (static_cast<int>(mpart.size()) != space_->natom()) {
      updateBase(mpart, flag, uptype, neigh_, neighOne_, neighOneOld_);
    }
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

void PairSquareWell::initHardSphere(const int itype, const int jtype) {
  rCutijset(itype, jtype, sigij_[itype][jtype]);
}

double PairSquareWell::allPartEnerForce(const int flag) {
  peSRone_ = 0;
  if (flag == 0) {
    peSRone_ = peTot();
    return peSRone_;
  } else {
    if ( (space_->cellType() == 0) || (rCutMaxAll_ > space_->dCellMin()) ) {
      return allPartEnerForceNoCell();
    } else if (space_->cellType() == 1) {
      return allPartEnerForceCell();
    } else {
      ASSERT(0, "cell type(" << space_->cellType() << ")");
    }
  }
  return 1e300;
}

double PairSquareWell::allPartEnerForceNoCell() {
  // shorthand for read-only space variables
  const int natom = space_->natom();
  const vector<double> x = space_->x();
  const vector<int> mol = space_->mol();
  const vector<int> type = space_->type();

  // declare variables for optimization
  double r2, xi, yi, zi;
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
        if ( (iMol != jMol) && (eps_[jtype] != 0) ) {
          // separation distance with periodic boundary conditions
          vector<double> xij(dimen_);
          xij[0] = xi - x[dimen_*jpart];
          xij[1] = yi - x[dimen_*jpart+1];
          xij[2] = zi - x[dimen_*jpart+2];
          vector<double> dxpbc = space_->pbc(xij);
          for (int dim = 0; dim < dimen_; ++dim) {
            xij[dim] += dxpbc[dim];
          }
          r2 = vecDotProd(xij, xij);

          // no interaction beyond cut-off distance
          const double rCutij = rCutij_[itype][jtype];
          if (r2 < rCutij*rCutij) {
            const double sigij = sigij_[itype][jtype];
            if (r2 < sigij*sigij) {
              // hard sphere if less than sigma
              peSRone_ += NUM_INF;
            } else {
              // otherwise, square well
              peSRone_ -= epsij_[itype][jtype];
            }
          }
        }
      }
    }
  }
  return peSRone_;
}

double PairSquareWell::allPartEnerForceCell() {
  // shorthand for read-only space variables
  const vector<double> &x = space_->x();
  const vector<int> &mol = space_->mol();
  const vector<int> &type = space_->type();
  const vector<vector<int> > neighCell = space_->neighCell();
  const vector<vector<int> > cellList = space_->cellList();

  // declare variables for optimization
  double r2, xi, yi, zi;
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
              if ( (iMol != jMol) && (eps_[jtype] != 0) &&
                   ((jCell != iCell) || (ipart < jpart)) ) {
                // separation distance with periodic boundary conditions
                vector<double> xij(dimen_);
                xij[0] = xi - x[dimen_*jpart];
                xij[1] = yi - x[dimen_*jpart+1];
                xij[2] = zi - x[dimen_*jpart+2];
                vector<double> dxpbc = space_->pbc(xij);
                for (int dim = 0; dim < dimen_; ++dim) {
                  xij[dim] += dxpbc[dim];
                }
                r2 = vecDotProd(xij, xij);

                // no interaction beyond cut-off distance
                const double rCutij = rCutij_[itype][jtype];
                if (r2 < rCutij*rCutij) {
                  const double sigij = sigij_[itype][jtype];
                  if (r2 < sigij*sigij) {
                    // hard sphere if less than sigma
                    peSRone_ += NUM_INF;
                  } else {
                    // otherwise, square well
                    peSRone_ -= epsij_[itype][jtype];
                  }
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

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

