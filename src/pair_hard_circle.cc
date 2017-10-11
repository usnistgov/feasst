/**
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of this agreement (see LICENSE.txt) and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#include <string>
#include "pair_hard_circle.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

PairHardCircle::PairHardCircle(Space* space,
  const double rCut)
  : Pair(space, rCut) {
  defaultConstruction_();
}
PairHardCircle::PairHardCircle(Space* space,
  const char* fileName)
  : Pair(space, fileName) {
  defaultConstruction_();
  dCircle_ = fstod("dCircle", fileName);
  rDep_ = fstod("rDep", fileName);
}

void PairHardCircle::defaultConstruction_() {
  className_.assign("PairHardCircle");
  dCircle_ = 1;
  initRDep();
  ASSERT(space_->dimen() == 2,
    "Hard circle potential requires spatial dimenion(" << space_->dimen()
    << ") of 2");
}

void PairHardCircle::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  file << "# dCircle " << dCircle_ << endl;
  file << "# rDep " << rDep_ << endl;
}

void PairHardCircle::initEnergy() {
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
  fill(0., f_);
  fill(0., vr_);
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
              peSR_ += NUM_INF;
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
}

double PairHardCircle::multiPartEner(
  const vector<int> mpart,
  const int flag) {
  peSRone_ = 0.;
  if (flag == 0) {}  // remove unused parameter warning
  return multiPartEnerAtomCut2D(mpart);
}

void PairHardCircle::update(
  const vector<int> mpart,
  const int flag,
  const char* uptype) {
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
          peSRone_ += NUM_INF;
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
    peSRone_ += NUM_INF;
  } else {
    const double r = sqrt(r2),
                 R = 0.5*dCircle_ + rDep_;
    peSRone_ -= (2*R*R*acos(r*0.5/R) - r*sqrt(R*R-r*r*0.25))/PI/rDep_/rDep_;
  }
}

shared_ptr<PairHardCircle> makePairHardCircle(Space* space, const double rCut) {
  return make_shared<PairHardCircle>(space, rCut);
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_


