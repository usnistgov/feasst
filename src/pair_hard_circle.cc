/**
 * \file
 *
 * \brief pairwise interactions between hard circles in depletant
 */

#include "pair_hard_circle.h"

namespace feasst {

/**
 * Constructor
 */
PairHardCircle::PairHardCircle(Space* space,
  const double rCut)  //!< interaction cut-off distance
  : Pair(space, rCut) {
  defaultConstruction();
}
PairHardCircle::PairHardCircle(Space* space,
  const char* fileName)
  : Pair(space, fileName) {
  defaultConstruction();
  dCircle_ = fstod("dCircle", fileName);
  rDep_ = fstod("rDep", fileName);
}

/**
 * defaults in constructor
 */
void PairHardCircle::defaultConstruction() {
  className_.assign("PairHardCircle");
  dCircle_ = 1;
  rDep_ = 0.;
  ASSERT(space_->dimen() == 2,
    "Hard circle potential requires spatial dimenion(" << space_->dimen()
    << ") of 2");
}

/**
 * write restart file
 */
void PairHardCircle::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  file << "# dCircle " << dCircle_ << endl;
  file << "# rDep " << rDep_ << endl;
}

/**
 * Lennard-Jones pair-wise force calculation
 */
int PairHardCircle::initEnergy() {
  // shorthand for read-only space variables
  const int dimen = space_->dimen();
  const int natom = space_->natom();
  const vector<double> l = space_->l();
  const vector<double> x = space_->x();
  const vector<int> &mol = space_->mol();
  const vector<int> &type = space_->type();

  double pePart;
  const double R = 0.5*dCircle_ + rDep_;

  // zero accumulators: potential energy, force, and virial
  std::fill(pe_.begin(), pe_.end(), 0.);
  myFill(0., f_);
  myFill(0., vr_);
  peSR_ = 0;

  // loop through nearest neighbor atom pairs
  for (int ipart = 0; ipart < natom - 1; ++ipart) {
    if (eps_[type[ipart]] != 0) {
      for (int jpart = ipart + 1; jpart < natom; ++jpart) {
        if ( (mol[ipart] != mol[jpart]) && (eps_[type[jpart]] != 0) ) {
          // separation vector, xij with periodic boundary conditions
          double r, r2 = 0;
          vector<double> xij(dimen);
          for (int i = 0; i < dimen; ++i) {
            r = x[dimen_*ipart+i] - x[dimen_*jpart+i];
            if (r >  0.5 * l[i]) r -= l[i];
            if (r < -0.5 * l[i]) r += l[i];
            xij[i] = r;
            r2 += r*r;
          }

          // no interaction beyond cut-off distance
          if (r2 < rCut_*rCut_) {
            // hard sphere
            if (r2 < dCircle_*dCircle_) {
              peSR_ += std::numeric_limits<double>::max()/1e10;
            } else {
              // energy
              r = sqrt(r2);
              pePart = -(2*R*R*acos(r*0.5/R)
                - r*sqrt(R*R-r*r*0.25))/PI/rDep_/rDep_;
              peSR_ += pePart;
              pe_[ipart] += pePart/2.;
              pe_[jpart] += pePart/2.;

              // force
              vector<double> fij(dimen);
              const double fPart = 0.;
              for (int i = 0; i < dimen; ++i) {
                fij[i] += fPart * xij[i];
                f_[ipart][i] += fij[i];
                f_[jpart][i] -= fij[i];

                // virial tensor
                for (int j = 0; j < dimen; ++j) {
                  vr_[ipart][i][j] += 0.5 * xij[j] * fij[i];
                  vr_[jpart][i][j] += 0.5 * xij[j] * fij[i];
                }
              }
            }
          }
        }
      }
    }
  }
  peTot_ = peSR_;
  return 0;
}

/**
 * Lennard-Jones potential energy contribution due to particles
 */
double PairHardCircle::multiPartEner(
  const vector<int> mpart,    //!< particles to calculate energy interactions
  const int flag) {     //!< place holder for other pair styles
  peSRone_ = 0.;
  if (flag == 0) {}  // remove unused parameter warning
  return multiPartEnerAtomCut2D(mpart);
}

/**
 * stores, restores or updates variables to avoid recompute of entire
 * configuration after every change
 */
void PairHardCircle::update(
  const vector<int> mpart,    //!< particles involved in move
  const int flag,         //!< type of move
  const char* uptype) {    //!< description of update type
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
      peSR_ += peSRone_ - deSR_;
      peTot_ += peSRone_ - deSR_;
    }
    if (flag == 2) {
      peSR_ -= deSR_;
      peTot_ -= deSR_;
    }
    if (flag == 3) {
      peSR_ += deSR_;
      peTot_ += deSR_;
    }
  }
}

/**
 * potential energy and forces of all particles
 *  if flag == 0, dummy calculation
 *  if flag == 1, all config calculation
 */
double PairHardCircle::allPartEnerForce(const int flag) {
  peSRone_ = 0.;
  if (flag == 0) {
    peSRone_ = peTot();
    return peSRone_;
  } else {
    return allPartEnerForceNoCell();
  }
  return 1e300;
}

/**
 * potential energy and forces of all particles
 */
double PairHardCircle::allPartEnerForceNoCell() {
  // shorthand for read-only space variables
  const int nMol = space_->nMol();
  const vector<double> &l = space_->l();
  const vector<double> &x = space_->x();
  const vector<int> &mol2part = space_->mol2part();
  const double lx = l[0], ly = l[1], halflx = lx/2., halfly = ly/2.;
  const double R = 0.5*dCircle_ + rDep_;
  double r, r2, dx, dy;

  // loop through pairs of molecules
  for (int iMol = 0; iMol < nMol - 1; ++iMol) {
    const int ipart = mol2part[iMol];
    const double xi = x[dimen_*ipart+0];
    const double yi = x[dimen_*ipart+1];
    for (int jMol = iMol + 1; jMol < nMol; ++jMol) {
      const int jpart = mol2part[jMol];

      // separation distance with periodic boundary conditions
      dx = xi - x[dimen_*jpart];
      dy = yi - x[dimen_*jpart+1];
      if (dx >  halflx) dx -= lx;
      if (dx < -halflx) dx += lx;
      if (dy >  halfly) dy -= ly;
      if (dy < -halfly) dy += ly;
      r2 = dx*dx + dy*dy;

      // no interaction beyond cut-off distance
      if (r2 < rCutSq_) {
        // hard sphere
        if (r2 < dCircle_*dCircle_) {
          peSRone_ += std::numeric_limits<double>::max()/1e10;
        } else {
          r = sqrt(r2);
          peSRone_ -= (2*R*R*acos(r*0.5/R)
            - r*sqrt(R*R-r*r*0.25))/PI/rDep_/rDep_;
        }
      }
    }
  }
  return peSRone_;
}

void PairHardCircle::multiPartEnerAtomCutInner(const double &r2,
  const int &itype, const int &jtype) {
  if (itype == jtype) {}  // remove unused parameter warning
  // hard sphere
  if (r2 < dCircle_*dCircle_) {
    peSRone_ += std::numeric_limits<double>::max()/1e10;
  } else {
    const double r = sqrt(r2),
                 R = 0.5*dCircle_ + rDep_;
    peSRone_ -= (2*R*R*acos(r*0.5/R) - r*sqrt(R*R-r*r*0.25))/PI/rDep_/rDep_;
  }
}

}  // namespace feasst


