/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./pair_lj_coul_ewald.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <numeric>
#include "./pair.h"
#include "./space.h"
#include "./functions.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

PairLJCoulEwald::PairLJCoulEwald(Space* space,
  const double rCut)
  : Pair(space, rCut) {
  defaultConstruction_();
}

PairLJCoulEwald::PairLJCoulEwald(Space* space,
  const char* fileName)
  : Pair(space, fileName) {
  defaultConstruction_();
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

void PairLJCoulEwald::defaultConstruction_() {
  className_.assign("PairLJCoulEwald");
  alpha = 5.6 / space_->minl();
  initAtomCut(0);
  skipEPS0_ = 0;
}

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

void PairLJCoulEwald::initEnergy() {
  peSRone_ = 0;
  peLJone_ = 0;
  peQRealone_ = 0;
  allPartEnerForce(2);
  peLJ_ = peLJone_;
  peQReal_ = peQRealone_;

  // reciprical (Fourier) space
  forcesFrr_();

  // standard long range corrections
  lrcConf_();

  selfAll_();
}

double PairLJCoulEwald::multiPartEnerReal(
  const vector<int> mpart,
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
  double r2inv, r6inv, r2, epsij, sigij, dx, dy, dz;
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
            //qi = q_[itype];
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
              multiPartEnerAtomCutInner(r2, itype, type[jsite]);
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
      peLRCone_ += lrcOne_(ipart);
    }
  }

  // fourier space contributions
  peQFrrone_ = 0;

  return peLJone_ + peLRCone_ + peQRealone_ + peQFrrone_;
}

void PairLJCoulEwald::k2maxset(
  const int k2max) {
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
          const double k2 = vecDotProd(kvect, kvect);
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
  selfAll_();
}

void PairLJCoulEwald::selfAll_() {
  vector<int> tmp;
  selfCorrect_(tmp);
  peQFrrSelf_ = peQFrrSelfone_;
}

void PairLJCoulEwald::selfCorrect_(
  vector<int> mpart) {
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
        double r = sqrt(vecDotProd(xij, xij));
        peQFrrSelfone_ += q_[type[ipart]]*q_[type[jpart]]*erf(alpha*r)/r;
      }
    }
  }
}

double PairLJCoulEwald::peTot() {
  return peLJ_ + peLRC_ + peQReal_ + peQFrr_ - peQFrrSelf_;
}

void PairLJCoulEwald::delPart(const int ipart) {
  eikrx_.erase(eikrx_.begin() + ipart);
  eikry_.erase(eikry_.begin() + ipart);
  eikrz_.erase(eikrz_.begin() + ipart);
  eikix_.erase(eikix_.begin() + ipart);
  eikiy_.erase(eikiy_.begin() + ipart);
  eikiz_.erase(eikiz_.begin() + ipart);
  ASSERT(0, "Depreciated function delPart(const int ipart)");
  delPartBase_(ipart);
}

void PairLJCoulEwald::delPart(
  const vector<int> mpart) {
  delPartBase_(mpart);
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
  addPartBase_();
}

void PairLJCoulEwald::forcesFrr_() {
  strucfacr_.resize(kexp_.size());
  strucfacrnew_.resize(strucfacr_.size());
  strucfaci_.resize(kexp_.size());
  strucfacinew_.resize(strucfaci_.size());
  multiPartEnerFrr(space_->listAtoms(), 1);
  update(space_->listAtoms(), 5, "update");
}

void PairLJCoulEwald::lrcConf_() {
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

double PairLJCoulEwald::lrcOne_(
  const int ipart) {
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

void PairLJCoulEwald::multiPartEnerFrr(
  const vector<int> mpart,
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

double PairLJCoulEwald::multiPartEner(const vector<int> multiPart,
  const int flag) {
  multiPartEnerReal(multiPart, flag);
  if (cheapEnergy_ == false) {
    multiPartEnerFrr(multiPart, flag);

    if (flag == 0 || flag == 1) {
      peQFrrSelfone_ = 0;
    } else if (flag == 2 || flag == 3) {
      selfCorrect_(multiPart);
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

void PairLJCoulEwald::update(
  const vector<int> mpart,
  const int flag,
  const char* uptype) {
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
  rCutijset(0, 0, rCut_);
}

void PairLJCoulEwald::initBulkSPCE(
  const double alphatmp,
  const int kmax) {
  initBulkSPCE();
  initKSpace(alphatmp, kmax);
}

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

void PairLJCoulEwald::initKSpace(const double alphatmp, const int k2max) {
  ASSERT(space_->minl() != 0,
    "box dimensions must be set before initializing kspace");
  alpha = alphatmp/space_->minl();
  k2maxset(k2max);
  erft_.init(alpha, rCut_);
  ASSERT(epsij().size() > 0, "Must Pair::initData() before Pair::initKSpace");
  rCutijset(0, 0, rCut_);
  initEnergy();
}

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

int PairLJCoulEwald::multiPartEnerAtomCutInner(const double &r2,
  const int &itype, const int &jtype) {
  const double epsij = epsij_[itype][jtype];
  double enlj;
  if (epsij != 0) {
    // energy and force prefactor
    const double sigij = sigij_[itype][jtype];
    const double r2inv = sigij*sigij/r2;
    const double r6inv = r2inv*r2inv*r2inv;
    enlj = 4.*epsij*(r6inv*(r6inv - 1.));
    // flj = 48.*epsij*(r6inv*r2inv*(r6inv - 0.5));
  } else {
    enlj = 0.;  // flj = 0;
  }
  peLJone_ += enlj;

  // charge interactions
  const double enq = q_[itype]*q_[jtype]*erft_.eval(r2);
//  const double fq = enq + q_[itype]*q_[jtype]
//            *(2.*alpha*exp(-alpha*alpha*r2)/sqrt(PI));

  peQRealone_ += enq;
  peSRone_ += enq + enlj;
  return 1;
}

void PairLJCoulEwald::pairSiteSite_(const int &iSiteType, const int &jSiteType,
  double * energy, double * force, int * neighbor, const double &dx,
  const double &dy, const double &dz) {
  const double r2 = dx*dx + dy*dy + dz*dz;
  const double epsij = epsij_[iSiteType][jSiteType];
  *energy = 0.;
  *force = 0.;
  *neighbor = 1;
  if (epsij != 0) {
    // energy and force prefactor
    const double sigij = sigij_[iSiteType][jSiteType];
    const double r2inv = sigij*sigij/r2;
    const double r6inv = r2inv*r2inv*r2inv;
    const double enlj = 4.*epsij*(r6inv*(r6inv - 1.));
    *energy += enlj;
    peLJone_ += enlj;
    *force += 48.*epsij*(r6inv*r2inv*(r6inv - 0.5));
  }

  // charge interactions
  const double enq = q_[iSiteType]*q_[jSiteType]*erft_.eval(r2);
  *energy += enq;
  peQRealone_ += enq;
  peSRone_ += *energy;
  *force += q_[iSiteType]*q_[jSiteType]
            *(2.*alpha*exp(-alpha*alpha*r2)/sqrt(PI));
}

double PairLJCoulEwald::pairLoopSite_(
  const vector<int> &siteList,
  const int noCell) {
  peLJone_ = 0.;
  peSRone_ = 0.;
  peQRealone_ = 0.;
  return Pair::pairLoopSite_(siteList, noCell);
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_



