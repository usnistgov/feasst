/**
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of this agreement (see LICENSE.txt) and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#include "./pair_lj.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

PairLJ::PairLJ(Space* space,
  const double rCut)
  : Pair(space, rCut) {
  defaultConstruction_();
}

PairLJ::PairLJ(Space* space,
  const char* fileName)
  : Pair(space, fileName) {
  defaultConstruction_();

  string strtmp = fstos("linearShiftFlag", fileName);
  if (!strtmp.empty()) {
    linearShift(stoi(strtmp));
  } else {
    string strtmp = fstos("cutShiftFlag", fileName);
    if (!strtmp.empty()) cutShift(stoi(strtmp));
  }
}

void PairLJ::defaultConstruction_() {
  className_.assign("PairLJ");
  peLRConePreCalc_ = (8./3.)*PI*((1./3.)*pow(rCut_, -9) - pow(rCut_, -3));
  linearShift();
}

void PairLJ::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  if (linearShiftFlag_ != 0) {
    file << "# linearShiftFlag " << linearShiftFlag_ << endl;
  } else if (cutShiftFlag_ != 0) {
    file << "# cutShiftFlag " << cutShiftFlag_ << endl;
  }
}

void PairLJ::initEnergy() {
  // shorthand for read-only space variables
  const int dimen = space_->dimen();
  const int natom = space_->natom();
  const vector<double> l = space_->l();
  const vector<double> x = space_->x();
  const vector<int> mol = space_->mol();
  const vector<int> type = space_->type();

  // zero accumulators: potential energy, force, and virial
  std::fill(pe_.begin(), pe_.end(), 0.);
  fill(0., f_);
  fill(0., vr_);
  peLJ_ = 0;
  fCOM_.clear();
  fCOM_.resize(space_->nMol(), vector<double>(dimen_, 0.));

  double r2inv, r6inv, pePart;

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
            // energy
            r2inv = 1./r2;
            r6inv = r2inv*r2inv*r2inv;
            pePart = 4. * (r6inv*(r6inv - 1.)) + peShift_;
            if (linearShiftFlag_) pePart += peLinearShift_ * (sqrt(r2) - rCut_);
            peLJ_ += pePart;
            pe_[ipart] += pePart/2.;
            pe_[jpart] += pePart/2.;

            // force
            vector<double> fij(dimen);
            const double fPart = 48. * (r6inv*r2inv*(r6inv - 0.5));
            for (int i = 0; i < dimen; ++i) {
              fij[i] += fPart * xij[i];
              f_[ipart][i] += fij[i];
              f_[jpart][i] -= fij[i];
              fCOM_[mol[ipart]][i] += fij[i];
              fCOM_[mol[jpart]][i] -= fij[i];

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

  // standard long range corrections
  peLRC_ = 0;
  if (lrcFlag == 1) {
    const double enlrc = (8./3.)*PI*(space_->natom()
      /space_->vol())*((1./3.)*pow(rCut_, -9) - pow(rCut_, -3));
    const double vrlrc =(16./3.)*PI*(space_->natom()
      /space_->vol())*((2./3.)*pow(rCut_, -9) - pow(rCut_, -3));
    for (int ipart = 0; ipart < natom; ++ipart) {
      pe_[ipart] += enlrc;
      peLRC_ += enlrc;
      for (int i = 0; i < dimen; ++i) {
        vr_[ipart][i][i] += vrlrc;
      }
    }
  }

  peTot_ = peLJ_ + peLRC_;
}

double PairLJ::multiPartEner(
  const vector<int> mpart,
  const int flag) {
  if (flag == 0) {}  // remove unused parameter warning

  // zero potential energy contribution of particle ipart
  peSRone_ = 0;
  peLRCone_ = 0;


  // contribution of particles to standard long range corrections
  // (non-additive ~n^2)
  if (lrcFlag == 1) {
    peLRCone_ += static_cast<int>(mpart.size())*(2.*space_->natom() - 1.)
      /space_->vol()*peLRConePreCalc_;
  }

  if (neighOn_) {
    ASSERT(mpart.size() <= 1, "assumes molecules of size 1");
    return multiPartEnerNeigh(mpart) + peLRCone_;
  } else {
    if (mpart.size() == 1) {
      return multiPartEnerNoNeigh(mpart) + peLRCone_;
    } else {
      if (atomCut_ == 1) {
        if (static_cast<int>(epsij_.size()) == space_->nParticleTypes()) {
          return multiPartEnerAtomCutNoNeigh(mpart) + peLRCone_;
        }
      } else {
        ASSERT(0, "molecular cutoff with multiple particles not implemented");
      }
    }
  }
  return 1e100;
}

double PairLJ::multiPartEnerAtomCutNoNeigh(
  const vector<int> mpart
  ) {
  // shorthand for read-only space variables
  const int natom = space_->natom();
  const vector<int> type = space_->type();
  const vector<double> &x = space_->x();
  const vector<double> &l = space_->l();
  const vector<int> &mol = space_->mol();

  // declare variables for optimization
  double r6inv, r2, xi, yi, zi, dx, dy, dz;
  const double lx = l[0], ly = l[1], lz = l[2],
    halflx = lx/2., halfly = ly/2., halflz = lz/2.;

  // loop through all particles in mpart
  for (unsigned int ii = 0; ii < mpart.size(); ++ii) {
    const int ipart = mpart[ii];
    if (eps_[type[ipart]] != 0) {   // ignore non-interacting atoms
      const int iMol = mol[ipart];
      xi = x[dimen_*ipart];
      yi = x[dimen_*ipart+1];
      zi = x[dimen_*ipart+2];

      // loop through all particles interacting with ipart
      for (int jpart = 0; jpart < natom; ++jpart) {
        const int jtype = type[jpart];
        if (eps_[jtype] != 0) {
          if (iMol != mol[jpart]) {
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
            if (r2 < rCutSq_) {
              r6inv = 1./(r2*r2*r2);
              peSRone_ += 4. * (r6inv*(r6inv - 1.)) + peShift_;
              if (linearShiftFlag_) {
                peSRone_ += peLinearShift_ * (sqrt(r2) - rCut_);
              }
            }
          }
        }
      }
    }
  }
  return peSRone_;
}

double PairLJ::multiPartEnerNeigh(
  const vector<int> mpart
  ) {
  // shorthand for read-only space variables
  const int natom = space_->natom();
  const vector<double> &x = space_->x();
  const vector<double> &l = space_->l();

  // declare variables for optimization
  int ipart, jpart, impart, mpartsize = static_cast<int>(mpart.size());
  double r6inv, r2, xi, yi, zi, dx, dy, dz, lx = l[0], ly = l[1], lz = l[2];
  double halflx = lx/2., halfly = ly/2., halflz = lz/2.;

  if (neighOn_) {
    neighOne_.clear();
    neighOne_.resize(static_cast<int>(mpart.size()), vector<int>());
  }

  for (impart = 0; impart < mpartsize; ++impart) {
    ipart = mpart[impart];
    xi = x[dimen_*ipart+0];
    yi = x[dimen_*ipart+1];
    zi = x[dimen_*ipart+2];

    for (jpart = 0; jpart < natom; ++jpart) {
      if (jpart != ipart) {
        // separation distance with periodic boundary conditions
        dx = xi - x[dimen_*jpart+0];
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
        if (r2 < rCutSq_) {
          // energy
          r6inv = 1./(r2*r2*r2);
          peSRone_ += 4. * (r6inv*(r6inv - 1.)) + peShift_;
          if (linearShiftFlag_) peSRone_ += peLinearShift_ * (sqrt(r2) - rCut_);

          // store new neighbor list
          if (neighOn_) {
            if ( (r2 < neighAboveSq_) && (r2 > neighBelowSq_) ) {
              neighOne_[impart].push_back(jpart);
            }
          }
        }
      }
    }
  }
  return peSRone_;
}

double PairLJ::multiPartEnerNoNeigh(
  const vector<int> mpart
  ) {
  // shorthand for read-only space variables
  const int natom = space_->natom();
  const vector<double> &x = space_->x();
  const vector<double> &l = space_->l();

  // declare variables for optimization
  int ipart, jpart, impart, mpartsize = static_cast<int>(mpart.size());
  double r6inv, r2, xi, yi, zi, dx, dy, dz, lx = l[0], ly = l[1], lz = l[2];
  double halflx = lx/2., halfly = ly/2., halflz = lz/2.;

  for (impart = 0; impart < mpartsize; ++impart) {
    ipart = mpart[impart];
    xi = x[dimen_*ipart+0];
    yi = x[dimen_*ipart+1];
    zi = x[dimen_*ipart+2];

    for (jpart = 0; jpart < natom; ++jpart) {
      if (jpart != ipart) {
        // separation distance with periodic boundary conditions
        dx = xi - x[dimen_*jpart+0];
        dy = yi - x[dimen_*jpart+1];
        dz = zi - x[dimen_*jpart+2];
//        slower version of periodic boundary conditions from arbitrary number
//        of box length wraps
//        dx += lx*rstatic_cast<int>(dx/lx); //
//        dy += ly*rstatic_cast<int>(dy/ly);
//        dz += lz*rstatic_cast<int>(dz/lz);
        if (dx >  halflx) dx -= lx;
        if (dx < -halflx) dx += lx;
        if (dy >  halfly) dy -= ly;
        if (dy < -halfly) dy += ly;
        if (dz >  halflz) dz -= lz;
        if (dz < -halflz) dz += lz;
        r2 = dx*dx + dy*dy + dz*dz;

        // no interaction beyond cut-off distance
        if (r2 < rCutSq_) {
          // energy
          r6inv = 1./(r2*r2*r2);
          peSRone_ += 4. * (r6inv*(r6inv - 1.)) + peShift_;
          if (linearShiftFlag_) peSRone_ += peLinearShift_ * (sqrt(r2) - rCut_);
        }
      }
    }
  }
  return peSRone_;
}

void PairLJ::update(const vector<int> mpart,
                    const int flag,
                    const char* uptype
  ) {
  if (neighOn_) {
    // rebuilt neighlist if all particles are updated
    if (static_cast<int>(mpart.size()) != space_->natom()) {
      updateBase(mpart, flag, uptype, neigh_, neighOne_, neighOneOld_);
  //  updateBase(mpart, flag, uptype, neighCut_, neighCutOne_, neighCutOneOld_);
    }
  }
  std::string uptypestr(uptype);

  if (uptypestr.compare("store") == 0) {
    if (flag == 0 || flag == 2 || flag == 3) {
      deLJ_ = peSRone_;
      deLRC_ = peLRCone_;
//      cout << "storing deLJ_ " << deLJ_ << " deLRC_ " << deLRC_ << endl;
    }
  }

  if (uptypestr.compare("update") == 0) {
    if (flag == 0) {
      peLJ_ += peSRone_ - deLJ_;
      peLRC_ += peLRCone_ - deLRC_;
      peTot_ += peSRone_ - deLJ_ + peLRCone_ - deLRC_;
//      cout << "de " << peSRone_ - deLJ_ << endl;
//      cout << "delrc " << peLRCone_ - deLRC_ << endl;
    }
    if (flag == 2) {
      peLJ_ -= deLJ_;
      peLRC_ -= deLRC_;
      peTot_ -= deLJ_ + deLRC_;
    }
    if (flag == 3) {
      peLJ_ += deLJ_;
      peLRC_ += deLRC_;
      peTot_ += deLJ_ + deLRC_;
    }
  }
  // cout << "update f" << flag << " pet " << peTot_ << endl;
}

void PairLJ::cutShift(const int flag
  ) {
  cutShiftFlag_ = flag;
  if (flag == 0) {
    peShift_ = 0;
  } else if (flag == 1) {
    peShift_ = -4.*(pow(rCut_, -12) - pow(rCut_, -6));
    lrcFlag = 0;
  } else {
    ASSERT(0, "Unrecognized cutShift flag(" << flag << ").");
  }
}

void PairLJ::linearShift(const int flag
  ) {
  cutShift(flag);
  linearShiftFlag_ = flag;
  if (flag == 0) {
    peLinearShift_ = 0;
  } else if (flag == 1) {
    peLinearShift_ = -4.*(-12.*pow(rCut_, -13) + 6.*pow(rCut_, -7));
  } else {
    ASSERT(0, "Unrecognized linearShift flag(" << flag << ").");
  }
}

shared_ptr<PairLJ> makePairLJ(Space* space, const double rCut) {
  return make_shared<PairLJ>(space, rCut);
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

