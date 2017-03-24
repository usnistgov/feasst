/**
 * \file
 *
 * \brief lennard-jones pairwise interactions with long range coulombic interactions treated by ewald
 *
 */

#include "./pair_lj_coul_ewald.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <numeric>
#include "./pair.h"
#include "./space.h"
#include "./functions.h"

/**
 * Constructor
 */
PairLJCoulEwald::PairLJCoulEwald(Space* space,
  const double rCut)  //!< interaction cut-off distance
  : Pair(space, rCut) {
  className_.assign("PairLJCoulEwald");
  alpha = 5.6 / space_->minl();
}
PairLJCoulEwald::PairLJCoulEwald(Space* space,
  const char* fileName)
  : Pair(space, fileName) {
  className_.assign("PairLJCoulEwald");
  k2max_ = fstoi("k2max", fileName);
  alpha  = fstod("alpha", fileName);
  const int nCharges = fstoi("nCharges", fileName);
  q_.resize(nCharges);
  for (int i = 0; i < nCharges; ++i) {
    stringstream ss;
    ss << "qCharge" << i;
    q_[i] = fstod(ss.str().c_str(), fileName);
  }
  initKSpace(alpha, k2max_);
}

PairLJCoulEwald::~PairLJCoulEwald() {
  sizeCheck();
  if (neighOn_) checkNeigh();
}

/**
 * write restart file
 */
void PairLJCoulEwald::writeRestart(const char* fileName) {
  Pair::writeRestart(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  file << "# k2max " << k2max_ << endl;
  file << "# alpha " << alpha*space_->minl() << endl;
  file << "# nCharges " << q_.size() << endl;
  for (unsigned int i = 0; i < q_.size(); ++i) {
    file << std::setprecision(std::numeric_limits<double>::digits10+2)
         << "# qCharge" << i << " " << q_[i] << endl;
  }
}

/**
 * Lennard-Jones pair-wise force calculation
 */
int PairLJCoulEwald::initEnergy() {
  // shorthand for read-only space variables
  const int dimen = space_->dimen();
  const int nMol = space_->nMol();
  const vector<double> x = space_->x();
  const vector<int> type = space_->type();
  const vector<int> mol2part = space_->mol2part();

  double r2inv, r6inv, enlj, flj;
  double enq, fq;
  peLJ_ = 0; peQReal_ = 0;

  // zero accumulators: potential energy, force, and virial
  std::fill(pe_.begin(), pe_.end(), 0.);
  myFill(0., f_);
  myFill(0., vr_);

  // loop through pairs of molecules
  for (int iMol = 0; iMol < nMol - 1; ++iMol) {
    const int ipart = mol2part[iMol];
    for (int jMol = iMol + 1; jMol < nMol; ++jMol) {
      const int jpart = mol2part[jMol];

      // separation vector, xij with periodic boundary conditions
      vector<double> xij(dimen);
      for (int i = 0; i < dimen; ++i) {
        xij[i] = x[dimen_*ipart+i] - x[dimen_*jpart+i];
      }
      const vector<double> dx = space_->pbc(xij);
      for (int i = 0; i < dimen; ++i) {
        xij[i] += dx[i];
      }

      // separation distance
      double r2 = 0.;
      for (vector<double>::iterator ptr = xij.begin();
           ptr != xij.end(); ++ptr) {
        r2 += (*ptr)*(*ptr);
      }

      // no interaction beyond molecular cut-off distance
      if (r2 < rCutSq_) {
        // loop through pairs of sites
        for (int isite = ipart; isite < mol2part[iMol+1]; ++isite) {
          for (int jsite = jpart; jsite < mol2part[jMol+1]; ++jsite) {
            // separation vector, xij with periodic boundary conditions
            vector<double> xij(dimen);
            for (int i = 0; i < dimen; ++i) {
              xij[i] = x[dimen_*isite+i] - x[dimen_*jsite+i];
            }
            const vector<double> dx = space_->pbc(xij);
            for (int i = 0; i < dimen; ++i) {
              xij[i] += dx[i];
            }

            // separation distance
            double r2 = 0.;
            for (vector<double>::iterator ptr = xij.begin();
                 ptr != xij.end(); ++ptr) {
              r2 += (*ptr)*(*ptr);
            }

            // dispersion interactions
            const double epsij = epsij_[type[isite]][type[jsite]];
            if ( epsij != 0 ) {
              // energy and force prefactor
              const double sigij = sigij_[type[isite]][type[jsite]];
              r2inv = sigij*sigij/r2;
              r6inv = r2inv*r2inv*r2inv;
              enlj = 4.*epsij*(r6inv*(r6inv - 1.));
              flj = 48.*epsij*(r6inv*r2inv*(r6inv - 0.5));
            } else {
              enlj = 0.; flj = 0;
            }
            peLJ_ += enlj;

            // charge interactions
            enq = q_[type[isite]]*q_[type[jsite]]*erft_.eval(r2);
//            r = sqrt(r2);
//            enq = q_[type[isite]]*q_[type[jsite]]*erfc(alpha*r)/r;
            fq = enq + q_[type[isite]]*q_[type[jsite]]
                      *(2.*alpha*exp(-alpha*alpha*r2)/sqrt(PI));

            peQReal_ += enq;
            pe_[isite] += (enlj + enq)/2.;
            pe_[jsite] += (enlj + enq)/2.;

            // force
            vector<double> fij(dimen);
            for (int i = 0; i < dimen; ++i) {
              fij[i] += (flj + fq) * xij[i];
              f_[isite][i] += fij[i];
              f_[jsite][i] -= fij[i];

              // virial tensor
              for (int j = 0; j < dimen; ++j) {
                vr_[isite][i][j] += 0.5 * xij[j] * fij[i];
                vr_[jsite][i][j] += 0.5 * xij[j] * fij[i];
              }
            }
          }
        }
      }
    }
  }

  // reciprical (Fourier) space
  forcesFrr();

  // standard long range corrections
  lrcConf();

  return 0;
}

/**
 * Real space energy contribution due to a subset of particles
 */
double PairLJCoulEwald::multiPartEnerReal(
  const vector<int> mpart,  //!< particles to calculate energy interactions
  const int flag) {
  // shorthand for read-only space variables
  const vector<double> &x = space_->x();
  const vector<double> l = space_->l();
  const vector<int> type = space_->type();
  const vector<int> mol2part = space_->mol2part();

  if (flag == 0) {}  // remove unused parameter warning

  // zero
  peLJone_ = 0;
  peLRCone_ = 0;
  peQRealone_ = 0;
  if (neighOn_) {
    neighOne_.clear();
    neighOne_.resize(1);
  }

  // optimization precomputations
  double r2inv, r6inv, r2, epsij, sigij, dx, dy, dz, qi, qij;
  const double lx = l[0], ly = l[1], lz = l[2],
    halflx = lx/2., halfly = ly/2., halflz = lz/2.;
  int jtype, itype, isite, jsite;
  const int ipart = mpart[0];
  const int iMol = space_->mol()[ipart];
  const double xi = x[dimen_*ipart+0];
  const double yi = x[dimen_*ipart+1];
  const double zi = x[dimen_*ipart+2];

  // use full list of atoms
  const vector<int> &neigh = space_->listMols();
  ASSERT(static_cast<int>(neigh.size()) == space_->nMol(),
    "problem with list of neighbors");

  // loop through molecules
  for (unsigned int ineigh = 0; ineigh < neigh.size(); ++ineigh) {
    const int jMol = neigh[ineigh];
    const int jpart = mol2part[jMol];
    if (iMol != jMol) {
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
        if (cheapEnergy_) {
          // dispersion interactions
          itype = type[ipart];
          jtype = type[jpart];
          epsij = epsij_[itype][jtype];
          sigij = sigij_[itype][jtype];
          r2inv = sigij*sigij/r2;
          r6inv = r2inv*r2inv*r2inv;
          peLJone_ += 4.*epsij*(r6inv*(r6inv - 1.));
        } else {
          // store new neighbor list
          if (neighOn_) {
            if ( (r2 < neighAboveSq_) && (r2 > neighBelowSq_) ) {
              neighOne_[0].push_back(jMol);
            }
          }
          // loop through pairs of sites
          for (isite = ipart; isite < mol2part[iMol+1]; ++isite) {
            itype = type[isite];
            qi = q_[itype];
            for (jsite = jpart; jsite < mol2part[jMol+1]; ++jsite) {
              // separation distance with periodic boundary conditions
              dx = x[dimen_*isite] - x[dimen_*jsite];
              dy = x[dimen_*isite+1] - x[dimen_*jsite+1];
              dz = x[dimen_*isite+2] - x[dimen_*jsite+2];
              if (dx >  halflx) dx -= lx;
              if (dx < -halflx) dx += lx;
              if (dy >  halfly) dy -= ly;
              if (dy < -halfly) dy += ly;
              if (dz >  halflz) dz -= lz;
              if (dz < -halflz) dz += lz;
              r2 = dx*dx + dy*dy + dz*dz;

              // dispersion interactions
              jtype = type[jsite];
              epsij = epsij_[itype][jtype];
              if ( epsij != 0 ) {
                sigij = sigij_[itype][jtype];
                r2inv = sigij*sigij/r2;
                r6inv = r2inv*r2inv*r2inv;
                peLJone_ += 4.*epsij*(r6inv*(r6inv - 1.));
              }

              // charge interactions
              qij = qi*q_[jtype];
              peQRealone_ += qij*erft_.eval(r2);
            }
          }
        }
      }
    }
  }

  // standard long range correction contribution
  if (cheapEnergy_ == false) {
    for (unsigned int impart = 0; impart < mpart.size(); ++impart) {
      const int ipart = mpart[impart];
      peLRCone_ += lrcOne(ipart);
    }
  }

  // fourier space contributions
  peQFrrone_ = 0;

  return peLJone_ + peLRCone_ + peQRealone_ + peQFrrone_;
}

/**
 * Real space energy contribution due to a subset of particles
 */
double PairLJCoulEwald::multiPartEnerRealAtomCut(
  const vector<int> mpart,  //!< particles to calculate energy interactions
  const int flag) {
  // shorthand for read-only space variables
  const int dimen = space_->dimen();
  const int natom = space_->natom();
  const vector<double> &x = space_->x();
  const vector<double> l = space_->l();
  const vector<int> type = space_->type();
  const vector<int> mol = space_->mol();

  if (flag == 0) {}  // remove unused parameter warning

  // zero
  double r, r2inv, r6inv, enlj;
  double enq, r2, epsij, sigij;
  int jtype;
  const double rCut2 = rCut_*rCut_;
  peLJone_ = 0;
  peLRCone_ = 0;
  peQRealone_ = 0;
  neighOne_.clear();
  neighOne_.resize(mpart.size(), vector<int>());
  myFill(0, neighOne_);
//  neighCutOne_.clear();
//  neighCutOne_.resize(mpart.size(), vector<int>());
//  myFill(0, neighCutOne_);

  for (unsigned int impart = 0; impart < mpart.size(); ++impart) {
    const int ipart = mpart[impart];
    const int itype = type[ipart], imol = mol[ipart];
    vector<double> xi;
    for (int dim = 0; dim < dimen; ++dim) xi.push_back(x[dimen_*ipart+dim]);
    const double qi = q_[itype];

//    // loop through nearest neighbor atom pairs
//    vector<int> neigh;
//    if ( (flag == 0) || (flag == 2) ) {
//
//    // use neighCut_ for up-to-date neighbors because there was no change in
//    // configuration
//    neigh = neighCut_[ipart];
//    } else {
//
//      // use full list of atoms
//      neigh = space_->listAtoms();
//      if (static_cast<int>(neigh.size()) != space_->natom()) {
//        std::cerr << "problemo68548689690567 " << std::endl;
//      }
//    }

    for (int jpart = 0; jpart < natom; ++jpart) {
//    for (int ineigh = 0; ineigh < static_cast<int>(neigh.size()); ++ineigh) {
//      const int jpart = neigh[ineigh];

      // no intramolecular interactions
      if (imol != mol[jpart]) {
        // separation distance with periodic boundary conditions
        r2 = 0.;

        for (int i = 0; i < dimen; ++i) {
          r = xi[i] - x[dimen_*jpart+i];
          if (r >  0.5 * l[i]) r -= l[i];
          if (r < -0.5 * l[i]) r += l[i];
          r2 += r*r;
        }

        // no interaction beyond cut-off distance
        if (r2 < rCut2) {
          // dispersion interactions
          jtype = type[jpart];
          epsij = epsij_[itype][jtype];
          if ( epsij != 0 ) {
            // energy
            sigij = sigij_[itype][jtype];
            r2inv = sigij*sigij/r2;
            r6inv = r2inv*r2inv*r2inv;
            enlj = 4.*epsij*(r6inv*(r6inv - 1.));
          } else {
            enlj = 0.;
          }
          peLJone_ += enlj;

          // charge interactions
          enq = qi*q_[jtype]*erft_.eval(r2);
          peQRealone_ += enq;

          // store new neighbor list
//          neighCutOne_[impart].push_back(jpart);
          if ( (r2 < neighAboveSq_) && (r2 > neighBelowSq_) ) {
            neighOne_[impart].push_back(jpart);
          }
        }
      }
    }
    // standard long range correction contribution of ipart
    peLRCone_ += lrcOne(ipart);
  }

  // fourier space contributions
  peQFrrone_ = 0;

  return peLJone_ + peLRCone_ + peQRealone_ + peQFrrone_;
}

/**
 * resize the wave vector arrays for a given maximum wave vector
 *  compute self interactions in fourier space
 */
void PairLJCoulEwald::k2maxset(
  const int k2max) {      //!< maximum wave vector magnitude cut off
  k2max_ = k2max;
  kmax_ = static_cast<int>(sqrt(k2max))+1;
  kxmax_ = kmax_ + 1;
  kymax_ = 2*kmax_ + 1;
  kzmax_ = 2*kmax_ + 1;
  eikrx_.resize(space_->natom()*kxmax_);
  eikry_.resize(space_->natom()*kymax_);
  eikrz_.resize(space_->natom()*kzmax_);
  eikix_.resize(space_->natom()*kxmax_);
  eikiy_.resize(space_->natom()*kymax_);
  eikiz_.resize(space_->natom()*kzmax_);

  // precompute wave vectors and prefactors
  kexp_.clear();
  k_.clear();
  vector<double> kvect(3);
  for (int kx = 0; kx <= kmax_; ++kx) {
    for (int ky = -kmax_; ky <= kmax_; ++ky) {
      for (int kz = -kmax_; kz <= kmax_; ++kz) {
        int k2 = kx*kx + ky*ky + kz*kz;
        if ( (k2 < k2max_) && (k2 != 0) ) {  // allen tildesley, srsw
        // if ( (k2 <= k2max_) && (k2 != 0) ) {  // gerhard
          kvect[0] = 2.*PI*kx/space_->l(0);
          kvect[1] = 2.*PI*ky/space_->l(1);
          kvect[2] = 2.*PI*kz/space_->l(2);
          const double k2 = myVecDotProd(kvect, kvect);
          double factor = 1.;
          if (kx != 0) factor = 2;
          kexp_.push_back(2.*PI*factor*exp(-k2/4./alpha/alpha)/k2
                          /space_->vol());
          k_.push_back(kx);
          k_.push_back(ky+kmax_);
          k_.push_back(kz+kmax_);
        }
      }
    }
  }

  // precompute self interaction energy
  vector<int> tmp;
  selfCorrect(tmp);
  peQFrrSelf_ = peQFrrSelfone_;
}

/**
 * computes self interaction energies to compensate for Ewald Sum
 *  of mpart particles
 *  if mpart is null, updates self interaction energy for entire system
 */
void PairLJCoulEwald::selfCorrect(
  vector<int> mpart) {   //!< compute self interactions of mpart
  // read only space variables
  const int dimen = space_->dimen();
  const int natom = space_->natom();
  const vector<double> x = space_->x();
  const vector<int> type = space_->type();
  const vector<int> mol = space_->mol();

  int npart = static_cast<int>(mpart.size());
  int ipart, jpart;
  peQFrrSelfone_ = 0;

  if (npart == 0) {
    for (int ipart = 0; ipart < natom; ++ipart) {
      mpart.push_back(ipart);
    }
    npart = static_cast<int>(mpart.size());
  }
  ASSERT(npart <= natom, "mpart is too large");

  // compute fourier space self interaction energy in Ewald sum
  for (int i = 0; i < npart; ++i) {
    ipart = mpart[i];
    peQFrrSelfone_ += alpha*q_[type[ipart]]*q_[type[ipart]]/sqrt(PI);
  }
  for (int i = 0; i < npart - 1; ++i) {
    ipart = mpart[i];
    for (int j = i + 1; j < npart; ++j) {
      jpart = mpart[j];
      if (mol[ipart] == mol[jpart]) {
        // separation vector, xij with periodic boundary conditions
        vector<double> xij(dimen);
        for (int dim = 0; dim < dimen; ++dim) {
          xij[dim] = x[dimen_*ipart+dim] - x[dimen_*jpart+dim];
        }
        const vector<double> dx = space_->pbc(xij);
        for (int i = 0; i < dimen; ++i) {
          xij[i] += dx[i];
        }

        // separation distance
        double r = sqrt(myVecDotProd(xij, xij));
        peQFrrSelfone_ += q_[type[ipart]]*q_[type[jpart]]*erf(alpha*r)/r;
      }
    }
  }
}

/**
 * Returns the total potential energy of the system
 */
double PairLJCoulEwald::peTot() {
  return peLJ_ + peLRC_ + peQReal_ + peQFrr_ - peQFrrSelf_;
}

/**
 * delete one particle
 */
void PairLJCoulEwald::delPart(const int ipart) {     //!< particle to delete
  eikrx_.erase(eikrx_.begin() + ipart);
  eikry_.erase(eikry_.begin() + ipart);
  eikrz_.erase(eikrz_.begin() + ipart);
  eikix_.erase(eikix_.begin() + ipart);
  eikiy_.erase(eikiy_.begin() + ipart);
  eikiz_.erase(eikiz_.begin() + ipart);
  mout_("warning", std::ostringstream().flush()
    << "Depreciated function delPart(const int ipart)");
  delPartBase(ipart);
}

/**
 * delete particles
 */
void PairLJCoulEwald::delPart(
  const vector<int> mpart) {     //!< particles to delete
  delPartBase(mpart);
  int natom = space_->natom();
  for (int i = mpart.size() - 1; i >= 0; --i) {
    int ipart = mpart[i];
    if (fastDel_) {
      const int jpart = space_->natom() - mpart.size() + i;
      for (int k = kxmax_ - 1; k >= 0; --k) {
        eikrx_[natom*k+ipart] = eikrx_[natom*k+jpart];
        eikix_[natom*k+ipart] = eikix_[natom*k+jpart];
      }
      for (int k = kymax_ - 1; k >= 0; --k) {
        eikry_[natom*k+ipart] = eikry_[natom*k+jpart];
        eikiy_[natom*k+ipart] = eikiy_[natom*k+jpart];
        eikrz_[natom*k+ipart] = eikrz_[natom*k+jpart];
        eikiz_[natom*k+ipart] = eikiz_[natom*k+jpart];
      }
      ipart = jpart;
    }
    for (int k = kxmax_ - 1; k >= 0; --k) {
      eikrx_.erase(eikrx_.begin() + natom*k+ipart);
      eikix_.erase(eikix_.begin() + natom*k+ipart);
    }
    for (int k = kymax_ - 1; k >= 0; --k) {
      eikry_.erase(eikry_.begin() + natom*k+ipart);
      eikiy_.erase(eikiy_.begin() + natom*k+ipart);
      eikrz_.erase(eikrz_.begin() + natom*k+ipart);
      eikiz_.erase(eikiz_.begin() + natom*k+ipart);
    }
    --natom;
  }
}

/**
 * add one particle
 */
void PairLJCoulEwald::addPart() {
  const int nAdd = space_->natom() - eikrx_.size()/kxmax_;
  for (int i = 0; i < nAdd; ++i) {
    int natomprev = eikrx_.size()/kxmax_;
    for (int k = kxmax_ - 1; k >= 0; --k) {
      eikrx_.insert(eikrx_.begin() + natomprev*k+natomprev, 0.);
      eikix_.insert(eikix_.begin() + natomprev*k+natomprev, 0.);
    }
    for (int k = kymax_ - 1; k >= 0; --k) {
      eikry_.insert(eikry_.begin() + natomprev*k+natomprev, 0.);
      eikiy_.insert(eikiy_.begin() + natomprev*k+natomprev, 0.);
      eikrz_.insert(eikrz_.begin() + natomprev*k+natomprev, 0.);
      eikiz_.insert(eikiz_.begin() + natomprev*k+natomprev, 0.);
    }
  }
  addPartBase();
}

/**
 * reciprical space forces, energy and virial
 *  fourier space sum over wave vectors for Ewald sum
 */
void PairLJCoulEwald::forcesFrr() {
  strucfacr_.resize(kexp_.size());
  strucfacrnew_.resize(strucfacr_.size());
  strucfaci_.resize(kexp_.size());
  strucfacinew_.resize(strucfaci_.size());
  multiPartEnerFrr(space_->listAtoms(), 1);
  update(space_->listAtoms(), 5, "update");
}

/**
 * compute standard long range contributions of all particles
 */
void PairLJCoulEwald::lrcConf() {
  // shorthand for read-only space variables
  const int natom = space_->natom();
  const int dimen = space_->dimen();
  const vector<int> type = space_->type();

  // standard long range corrections for oxygen atoms
  peLRC_ = 0;
  if (lrcFlag == 1) {
    vector<int> nLJ(eps_.size());
    for (int ipart = 0; ipart < natom; ++ipart) {
      ++nLJ[type[ipart]];
    }
    for (int ipart = 0; ipart < natom; ++ipart) {
      const int itype = type[ipart];
      if (sig_[itype] != 0) {
        const double enlrc = (8./3.)*PI*(nLJ[itype]/space_->vol())
          *eps_[itype]*pow(sig_[itype], 3)*((1./3.)*pow(rCut_/sig_[itype], -9)
          - pow(rCut_/sig_[itype], -3));
        const double vrlrc =(16./3.)*PI*(nLJ[itype]/space_->vol())
          *eps_[itype]*pow(sig_[itype], 3)*((2./3.)*pow(rCut_/sig_[itype], -9)
          - pow(rCut_/sig_[itype], -3));
        pe_[ipart] += enlrc;
        peLRC_ += enlrc;
        for (int i = 0; i < dimen; ++i) {
          vr_[ipart][i][i] += vrlrc;
        }
      }
    }
  }
}

/**
 * compute standard long range contributions of one particle ipart
 */
double PairLJCoulEwald::lrcOne(
  const int ipart) {    //!< contributions from ipart
  // shorthand for read-only space variables
  const vector<int> type = space_->type();
  const vector<int> nType = space_->nType();
  const int itype = type[ipart];
  double lrcone = 0;
  const double sigrinv3 = pow(rCut_/sig_[itype], -3),
                   sig3 = pow(sig_[itype], 3);

  // contribution of one particle to standard long range corrections
  // (non-additive ~n^2) if oxygen
  if (lrcFlag == 1) {
    if (sig_[itype] != 0) {
      lrcone = (2.*nType[itype] - 1.)*(8./3.)*PI*eps_[itype]*sig3/space_->vol()
               *sigrinv3*((1./3.)*sigrinv3*sigrinv3 - 1.);
    } else {
      lrcone = 0.;
    }
  }
  return lrcone;
}

/**
 * reciprical space energy of multiple particles
 *  flag==0, energy of old configuration (no moves, ins or dels)
 *  flag==1, new configuration, same number of particles
 *  flag==2, old configuration, preparing to delete (same as 0)
 *  flag==3, just inserted particle
 */
void PairLJCoulEwald::multiPartEnerFrr(
  const vector<int> mpart,    //!< particle contribution
  const int flag) {
  if (flag == 0) {
    peQFrrone_ = peQFrr_;
    return;
  }

  // shorthand for read-only space variables
  const int natom = space_->natom();
  const vector<double> x = space_->x();
  const vector<int> type = space_->type();
  const vector<double> l = space_->l();
  const double twopilxi = 2.*PI/l[0],
               twopilyi = 2.*PI/l[1],
               twopilzi = 2.*PI/l[2];

  const int msize = mpart.size();

  eikrxnew_.clear();
  eikrynew_.clear();
  eikrznew_.clear();
  eikrxnew_.resize(kxmax_*msize);
  eikrynew_.resize(kymax_*msize);
  eikrznew_.resize(kzmax_*msize);
  eikixnew_.clear();
  eikiynew_.clear();
  eikiznew_.clear();
  eikixnew_.resize(kxmax_*msize);
  eikiynew_.resize(kymax_*msize);
  eikiznew_.resize(kzmax_*msize);
  for (unsigned int k = 0; k < strucfacr_.size(); ++k) {
    strucfacrnew_[k]  = strucfacr_[k];
  }
  for (unsigned int k = 0; k < strucfaci_.size(); ++k) {
    strucfacinew_[k]  = strucfaci_[k];
  }

  // compute new structure factor components when inserting or moving particles
  if (flag == 1 || flag == 3) {
    // calculate eikr of kx = 0, 1 and -1 explicitly
    for (int i = 0; i < msize; ++i) {
      const int ipart = mpart[i];
      eikrxnew_[i] = 1.; eikixnew_[i] = 0.;
      eikrynew_[msize*kmax_+i] = 1.; eikiynew_[msize*kmax_+i] = 0.;
      eikrznew_[msize*kmax_+i] = 1.; eikiznew_[msize*kmax_+i] = 0.;
      eikrxnew_[msize+i] = cos(twopilxi*x[dimen_*ipart]);
      eikixnew_[msize+i] = sin(twopilxi*x[dimen_*ipart]);
      eikrynew_[msize*(1+kmax_)+i] = cos(twopilyi*x[dimen_*ipart+1]);
      eikiynew_[msize*(1+kmax_)+i] = sin(twopilyi*x[dimen_*ipart+1]);
      eikrznew_[msize*(1+kmax_)+i] = cos(twopilzi*x[dimen_*ipart+2]);
      eikiznew_[msize*(1+kmax_)+i] = sin(twopilzi*x[dimen_*ipart+2]);
      eikrynew_[msize*(-1+kmax_)+i] = eikrynew_[msize*(1+kmax_)+i];
      eikrznew_[msize*(-1+kmax_)+i] = eikrznew_[msize*(1+kmax_)+i];
      eikiynew_[msize*(-1+kmax_)+i] = -eikiynew_[msize*(1+kmax_)+i];
      eikiznew_[msize*(-1+kmax_)+i] = -eikiznew_[msize*(1+kmax_)+i];
    }

    // compute remaining eikr by recursion relation
    for (int kx = 2; kx <= kmax_; ++kx) {
      for (int i = 0; i < msize; ++i) {
        eikrxnew_[msize*kx+i] = eikrxnew_[msize*(kx -1)+i]*eikrxnew_[msize*1+i]
          - eikixnew_[msize*(kx -1)+i]*eikixnew_[msize*(1)+i];
        eikixnew_[msize*(kx)+i] = eikrxnew_[msize*(kx -1)+i]
          *eikixnew_[msize*(1)+i] + eikixnew_[msize*(kx -1)+i]
          *eikrxnew_[msize*1+i];
      }
    }
    for (int ky = 2; ky <= kmax_; ++ky) {
      for (int i = 0; i < msize; ++i) {
        eikrynew_[msize*(ky+kmax_)+i] = eikrynew_[msize*(ky+kmax_-1)+i]
        *eikrynew_[msize*(1+kmax_)+i] - eikiynew_[msize*(ky+kmax_-1)+i]
        *eikiynew_[msize*(1+kmax_)+i];
        eikiynew_[msize*(ky+kmax_)+i] = eikrynew_[msize*(ky+kmax_-1)+i]
        *eikiynew_[msize*(1+kmax_)+i] + eikiynew_[msize*(ky+kmax_-1)+i]
        *eikrynew_[msize*(1+kmax_)+i];
        eikrynew_[msize*(-ky+kmax_)+i] = eikrynew_[msize*(ky+kmax_)+i];
        eikiynew_[msize*(-ky+kmax_)+i] = -eikiynew_[msize*(ky+kmax_)+i];
      }
    }
    for (int kz = 2; kz <= kmax_; ++kz) {
      for (int i = 0; i < msize; ++i) {
        eikrznew_[msize*(kz+kmax_)+i] = eikrznew_[msize*(kz+kmax_-1)+i]
        *eikrznew_[msize*(1+kmax_)+i] - eikiznew_[msize*(kz+kmax_-1)+i]
        *eikiznew_[msize*(1+kmax_)+i];
        eikiznew_[msize*(kz+kmax_)+i] = eikrznew_[msize*(kz+kmax_-1)+i]
        *eikiznew_[msize*(1+kmax_)+i] + eikiznew_[msize*(kz+kmax_-1)+i]
        *eikrznew_[msize*(1+kmax_)+i];
        eikrznew_[msize*(-kz+kmax_)+i] = eikrznew_[msize*(kz+kmax_)+i];
        eikiznew_[msize*(-kz+kmax_)+i] = -eikiznew_[msize*(kz+kmax_)+i];
      }
    }
  }

  peQFrrone_ = 0.;
  double eikrrold, eikrrnew, eikriold, eikrinew;
  double eikrx, eikix, eikry, eikiy, eikrz, eikiz;
  double eikrxnew, eikixnew, eikrynew, eikiynew, eikrznew, eikiznew;
  int kx, ky, kz, kdim;
  for (unsigned int k = 0; k < strucfacr_.size(); ++k) {
    kdim = dimen_*k;
    kx = k_[kdim];
    ky = k_[kdim+1];
    kz = k_[kdim+2];

    // structure factor
    for (int i = 0; i < msize; ++i) {
      const int ipart = mpart[i];
      eikrx = eikrx_[natom*kx+ipart];
      eikix = eikix_[natom*kx+ipart];
      eikry = eikry_[natom*ky+ipart];
      eikiy = eikiy_[natom*ky+ipart];
      eikrz = eikrz_[natom*kz+ipart];
      eikiz = eikiz_[natom*kz+ipart];

      eikrrold = eikrx*eikry*eikrz
               - eikix*eikiy*eikrz
               - eikix*eikry*eikiz
               - eikrx*eikiy*eikiz;
      eikriold = -eikix*eikiy*eikiz
               + eikrx*eikry*eikiz
               + eikrx*eikiy*eikrz
               + eikix*eikry*eikrz;

      eikrxnew = eikrxnew_[msize*kx+i];
      eikixnew = eikixnew_[msize*kx+i];
      eikrynew = eikrynew_[msize*ky+i];
      eikiynew = eikiynew_[msize*ky+i];
      eikrznew = eikrznew_[msize*kz+i];
      eikiznew = eikiznew_[msize*kz+i];

      eikrrnew = eikrxnew*eikrynew*eikrznew
               - eikixnew*eikiynew*eikrznew
               - eikixnew*eikrynew*eikiznew
               - eikrxnew*eikiynew*eikiznew;
      eikrinew = -eikixnew*eikiynew*eikiznew
               + eikrxnew*eikrynew*eikiznew
               + eikrxnew*eikiynew*eikrznew
               + eikixnew*eikrynew*eikrznew;

      strucfacrnew_[k] += q_[type[ipart]]*(eikrrnew - eikrrold);
      strucfacinew_[k] += q_[type[ipart]]*(eikrinew - eikriold);
    }
    peQFrrone_ += kexp_[k]*(strucfacrnew_[k]*strucfacrnew_[k]
                          + strucfacinew_[k]*strucfacinew_[k]);
  }
}

/**
 *  compute interaction contributions of multiple particles
 *  flag==0, energy of old configuration (no moves, ins or dels)
 *  flag==1, new configuration, same number of particles
 *  flag==2, old configuration, preparing to delete (same as 0)
 *  flag==3, just inserted particle
 */
double PairLJCoulEwald::multiPartEner(const vector<int> multiPart,
  const int flag) {
  multiPartEnerReal(multiPart, flag);
  if (cheapEnergy_ == false) {
    multiPartEnerFrr(multiPart, flag);

    if (flag == 0 || flag == 1) {
      peQFrrSelfone_ = 0;
    } else if (flag == 2 || flag == 3) {
      selfCorrect(multiPart);
    }
  }

  double detot = 0;
  if (flag == 0 || flag == 1) {
    if (cheapEnergy_) {
      detot = peLJone_;
    } else {
      detot = peLJone_ + peLRCone_ + peQRealone_ + peQFrrone_ - peQFrrSelfone_;
    }
  } else if (flag == 2) {
    if (cheapEnergy_) {
      detot = peLJone_;
    } else {
      detot = peLJone_ + peLRCone_ + peQRealone_
            - peQFrrone_ + peQFrr_ - peQFrrSelfone_;
    }
  } else if (flag == 3) {
    if (cheapEnergy_) {
      detot = peLJone_;
    } else {
      detot = peLJone_ + peLRCone_ + peQRealone_
            + peQFrrone_ - peQFrr_ - peQFrrSelfone_;
    }
  }

  return detot;
}

/**
 * stores, restores or updates variables to avoid recompute of entire
 * configuration after every change
 *  flag==0, energy of old configuration (no moves, ins or dels)
 *  flag==1, new configuration, same number of particles
 *  flag==2, old configuration, preparing to delete (same as 0)
 *  flag==3, just inserted particle
 *  flag==5, computing entire configuration
 */
void PairLJCoulEwald::update(
  const vector<int> mpart,    //!< particles involved in move
  const int flag,         //!< type of move
  const char* uptype) {    //!< description of update type
  if (neighOn_) {
    updateBase(mpart, flag, uptype, neigh_, neighOne_, neighOneOld_);
//  updateBase(mpart, flag, uptype, neighCut_, neighCutOne_, neighCutOneOld_);
  }
  std::string uptypestr(uptype);

  if (uptypestr.compare("store") == 0) {
    if (flag == 0 || flag == 2 || flag == 3) {
      deLJ_ = peLJone_;
      deLRC_ = peLRCone_;
      deQReal_ = peQRealone_;
      deQFrr_ = peQFrrone_;
      deQFrrSelf_ = peQFrrSelfone_;
    }
  }

  if (uptypestr.compare("update") == 0) {
    peQFrr_ = peQFrrone_;
    if (flag == 0) {
      peLJ_ += peLJone_ - deLJ_;
      peLRC_ += peLRCone_ - deLRC_;
      peQReal_ += peQRealone_ - deQReal_;
      peQFrrSelf_ += peQFrrSelfone_ - deQFrrSelf_;
    }
    if (flag == 2) {
      peLJ_ -= deLJ_;
      peLRC_ -= deLRC_;
      peQReal_ -= deQReal_;
      peQFrrSelf_ -= deQFrrSelf_;
    }
    if (flag == 3) {
      peLJ_ += deLJ_;
      peLRC_ += deLRC_;
      peQReal_ += deQReal_;
      peQFrrSelf_ += deQFrrSelf_;
    }
    if (flag == 0 || flag == 1 || flag == 3 || flag == 5) {
      const int msize = mpart.size();
      const int natom = space_->natom();
      for (int i = 0; i < msize; ++i) {
        const int ipart = mpart[i];
        for (int k = 0; k < kxmax_; ++k) {
          eikrx_[natom*k+ipart] = eikrxnew_[msize*k+i];
          eikix_[natom*k+ipart] = eikixnew_[msize*k+i];
        }
        for (int k = 0; k < kymax_; ++k) {
          eikry_[natom*k+ipart] = eikrynew_[msize*k+i];
          eikrz_[natom*k+ipart] = eikrznew_[msize*k+i];
          eikiy_[natom*k+ipart] = eikiynew_[msize*k+i];
          eikiz_[natom*k+ipart] = eikiznew_[msize*k+i];
        }
      }
    }
    if (flag == 0 || flag == 1 || flag == 2 || flag == 3 || flag == 5) {
      for (unsigned int k = 0; k < strucfacr_.size(); ++k) {
        strucfacr_[k] = strucfacrnew_[k];
      }
      for (unsigned int k = 0; k < strucfaci_.size(); ++k) {
        strucfaci_[k] = strucfacinew_[k];
      }
    }
  }
}

/**
 * initialize bulk simulation of SPCE waters, atoms listed as oxygen then two
 * hydrogens
 */
void PairLJCoulEwald::initBulkSPCE() {
  vector<double> eps(2), sig(2);
  // permitivity of free space, e^2*mol/kJ/A
  const double permitivity = permitivityVacuum/elementaryCharge/elementaryCharge
                            *1e3/1e10/avogadroConstant;
  const double qh = 0.4238/sqrt(4. * PI * permitivity);
  q_.resize(2);
  q_[1] = qh; q_[0] = -2.*qh;
  eps[1] = 0.; eps[0] = 0.650169581;
  sig[1] = 0.; sig[0] = 3.16555789;
  initPairParam(eps, sig);
  initAtomCut(0);
}

/**
 * initialize bulk simulation of SPCE waters, atoms listed as oxygen then two
 * hydrogens
 *  also set alpha and kappa parameters
 */
void PairLJCoulEwald::initBulkSPCE(
  const double alphatmp,    //!< Ewald alpha parameter
  const int kmax) {         //!< Ewald kmax
  initBulkSPCE();
  initKSpace(alphatmp, kmax);
}

/**
 * check size of class variables
 */
void PairLJCoulEwald::sizeCheck() {
  bool er = false;
  std::ostringstream ermesg;
  if (strucfacr_.size() != strucfacrnew_.size()) {
    er = true;
    ermesg << "strucfacr_ (" << strucfacr_.size() << ") != strucfacrnew_ ("
           << strucfacrnew_.size() << ") ";
  }
  if (strucfacr_.size() != kexp_.size()) {
    er = true;
    ermesg << "strucfacr_ (" << strucfacr_.size() << ") != kexp_ ("
           << kexp_.size() << ") ";
  }
  if (er) {
    ASSERT(0, "size check failure" << endl << ermesg.str());
  }
}

/**
 * tabular error function constructor
 */
erftable::erftable() {
}

/**
 * initialize table for f(x) = erfc(alpha*r)/r
 */
void erftable::init(const double alpha, const double rCut) {
  n_ = 2e5;
  ds_ = pow(2*rCut, 2)/static_cast<double>(n_);
  for (int i = 0; i < n_; ++i) {
    const double s = static_cast<double>(i) * ds_;
    const double rij = sqrt(s);
    vtab_.push_back(erfc(alpha*rij)/rij);
  }
}

/**
 * evaluate tabular eror function
 */
double erftable::eval(const double x) const {
  const double sds = x / ds_;
  const int k = static_cast<int>(sds);
  const double xi = sds - static_cast<double>(k);
  const double vk = vtab_[k], vk1 = vtab_[k+1], vk2 = vtab_[k+2];
  const double t1 = vk + (vk1 - vk) * xi;
  const double t2 = vk1 + (vk2 - vk1) * (xi - 1.);
  return t1 + (t2 - t1)*xi*0.5;
}


/**
 * initialize kspace, number of wave vectors and screening parameters
 */
void PairLJCoulEwald::initKSpace(const double alphatmp, const int k2max) {
  ASSERT(space_->minl() != 0,
    "box dimensions must be set before initializing kspace");
  alpha = alphatmp/space_->minl();
  k2maxset(k2max);
  erft_.init(alpha, rCut_);
  initEnergy();
}

/***
 * initialize with LAMMPS data file
 */
void PairLJCoulEwald::initLMPData(const string fileName) {
  Pair::initLMPData(fileName);

  // open LAMMPS data file
  std::ifstream file(fileName.c_str());
  ASSERT(file.good(), "cannot find lammps DATA file " << fileName);
  ASSERT(dimen_ == 3, "dimenions(" << dimen_ << ") must be 3");

  shared_ptr<Space> s = space_->findAddMolInList(fileName);

  // read until Atoms
  readUntil("Atoms", file);

  // read charges
  q_.resize(s->nParticleTypes());
  const double permitivity = permitivityVacuum/elementaryCharge/elementaryCharge
    *1e3/1e10/avogadroConstant;  // permitivity of free space, e^2*mol/kJ/A
  // assumes that epsilon is in units of kJ/mol
  for (int iAtom = 0; iAtom < s->natom(); ++iAtom) {
    int tmp, type;
    double qtmp;
    string line;
    file >> tmp >> tmp >> type >> qtmp;
    type -= 1;
    q_[type] = qtmp/sqrt(4.*PI*permitivity);
    std::getline(file, line);
  }
}



