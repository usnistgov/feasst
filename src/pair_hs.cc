/**
 * \file
 *
 * \brief hard sphere pairwise interactions
 *
 * Implementation of pair_hs class
 */

#include "./pair_hs.h"

/**
 * Constructor
 */
PairHS::PairHS(Space* space,
         const double rCut  //!< interaction cut-off distance
  )
    : Pair(space, rCut) {
  defaultConstruction();
}
PairHS::PairHS(Space* space,
         const char* fileName
  )
    : Pair(space, fileName) {
  defaultConstruction();
}

/**
 * defaults in constructor
 */
void PairHS::defaultConstruction() {
  className_.assign("PairHS");
}

/**
 * write restart file
 */
void PairHS::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
}

/**
 * pair-wise force calculation
 */
int PairHS::initEnergy() {
//  cout << "in forces" << endl;

  // shorthand for read-only space variables
  const int dimen = space_->dimen();
  const int natom = space_->natom();
  const vector<double> l = space_->l();
  const vector<double> x = space_->x();
  const vector<int> mol = space_->mol();
  const vector<int> type = space_->type();

  // zero accumulators: potential energy, force, and virial
  peTot_ = 0;

  // loop through nearest neighbor atom pairs
  for (int ipart = 0; ipart < natom - 1; ++ipart) {
    if (eps_[type[ipart]] != 0) {
      for (int jpart = ipart + 1; jpart < natom; ++jpart) {
        if ( intraCheck_(ipart, jpart, mol[ipart], mol[jpart]) &&
             (eps_[type[jpart]] != 0) ) {
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

//          cout << "r2 " << r2 << " i " << ipart << " j " << jpart << endl;

          // no interaction beyond cut-off distance
          double sigij = 0;
          if (sigij_.size() != 0) {
            sigij = sigij_[type[ipart]][type[jpart]];
          } else {
            sigij = rCut_;
          }
          if (r2 < sigij*sigij) {
          // if (r2 < rCut_*rCut_) {

            // energy
            peTot_ += 1e12;
          }
//          std::streamsize ss = cout.precision();
//          cout<< std::setprecision(std::numeric_limits<long double>::digits10+2) << "peSR " << peTot_ << " for ipart " << ipart << " jpart " << jpart << " eps " << epsij_[type[ipart]][type[jpart]] << " sig " << sigij_[type[ipart]][type[jpart]] << " r2 " << r2 << endl;
//          cout << std::setprecision(ss);
        }
      }
    }
  }

  return 0;
}

/**
 * potential energy contribution due to particles
 */
double PairHS::multiPartEner(
  const vector<int> mpart,    //!< particles to calculate energy interactions
  const int flag     //!< place holder for other pair styles
  ) {
  if (flag == 0) {}  // remove unused parameter warning

  // zero potential energy contribution of particle ipart
  peSRone_ = 0;
  return multiPartEnerAtomCut(mpart);
}

void PairHS::multiPartEnerAtomCutInner(const double &r2,
  const int &itype,
  const int &jtype) {
  const double sigij = sigij_[itype][jtype];
  if (r2 < sigij*sigij) {
    peSRone_ += 1e12;
  }
  // cout << "in Inner: r2 " << r2 << " itype " << itype << " jtype " << jtype
  //   << " peSRone " << peSRone_ << " sigij " << sigij << endl;
}

/**
 * stores, restores or updates variables to avoid recompute of entire
 * configuration after every change
 */
void PairHS::update(const vector<int> mpart,    //!< particles involved in move
                    const int flag,         //!< type of move
                    const char* uptype    //!< description of update type
  ) {
  if (neighOn_) {
    updateBase(mpart, flag, uptype, neigh_, neighOne_, neighOneOld_);
//  updateBase(mpart, flag, uptype, neighCut_, neighCutOne_, neighCutOneOld_);
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




