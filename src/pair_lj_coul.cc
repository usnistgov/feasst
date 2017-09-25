/**
 * \file
 *
 * \brief lennard-jones pairwise interactions with long range coulombic interactions without Ewald.
 *
 */

#include "./pair_lj_coul.h"
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

PairLJCoul::PairLJCoul(Space* space,
  const double rCut)  //!< interaction cut-off distance
  : Pair(space, rCut) {
  className_.assign("PairLJCoul");
}
PairLJCoul::PairLJCoul(Space* space,
  const char* fileName)
  : Pair(space, fileName) {
  className_.assign("PairLJCoul");
  const int nCharges = fstoi("nCharges", fileName);
  q_.resize(nCharges);
  for (int i = 0; i < nCharges; ++i) {
    stringstream ss;
    ss << "qCharge" << i;
    q_[i] = fstod(ss.str().c_str(), fileName);
  }
  string strtmp = fstos("hardSphereCOM", fileName);
  if (!strtmp.empty()) {
    hardSphereCOM = stod(strtmp);
  }
}

PairLJCoul::~PairLJCoul() {
  sizeCheck();
  if (neighOn_) checkNeigh();
}

/**
 * write restart file
 */
void PairLJCoul::writeRestart(const char* fileName) {
  Pair::writeRestart(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  file << "# nCharges " << q_.size() << endl;
  for (unsigned int i = 0; i < q_.size(); ++i) {
    file << std::setprecision(std::numeric_limits<double>::digits10+2)
         << "# qCharge" << i << " " << q_[i] << endl;
  }
  file << "# hardSphereCOM " << hardSphereCOM << endl;
}

/**
 * Lennard-Jones pair-wise force calculation
 */
void PairLJCoul::initEnergy() {
  // shorthand for read-only space variables
  const int dimen = space_->dimen();
  const int nMol = space_->nMol();
  const vector<double> x = space_->x();
  const vector<int> type = space_->type();
  const vector<int> mol2part = space_->mol2part();

  double r2inv, r6inv, enlj, flj;
  double enq;
  peLJ_ = 0; peQReal_ = 0;

  // zero accumulators: potential energy, force, and virial
  std::fill(pe_.begin(), pe_.end(), 0.);
  fill(0., f_);
  fill(0., vr_);

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
        // hard sphere interaction if below threshold
        if (r2 < hardSphereCOM*hardSphereCOM) {
          peLJ_ += NUM_INF;
        } else {
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
              enq = q_[type[isite]]*q_[type[jsite]]/sqrt(r2);

              peQReal_ += enq;
              pe_[isite] += (enlj + enq)/2.;
              pe_[jsite] += (enlj + enq)/2.;

              // force
              vector<double> fij(dimen);
              for (int i = 0; i < dimen; ++i) {
                fij[i] += flj * xij[i];
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
  }

  // standard long range corrections
  lrcConf();
}

/**
 * Real space energy contribution due to a subset of particles
 */
double PairLJCoul::multiPartEnerReal(
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
        if (r2 < hardSphereCOM*hardSphereCOM) {
          peLJone_ += NUM_INF;
        } else if (cheapEnergy_) {
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
              peQRealone_ += qij/sqrt(r2);
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

  return peLJone_ + peLRCone_ + peQRealone_;
}

/**
 * Real space energy contribution due to a subset of particles
 */
double PairLJCoul::multiPartEnerRealAtomCut(
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
  fill(0, neighOne_);
//  neighCutOne_.clear();
//  neighCutOne_.resize(mpart.size(), vector<int>());
//  fill(0, neighCutOne_);

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
          // hard sphere interaction if below threshold
          if (r2 < hardSphereCOM*hardSphereCOM) {
            peLJ_ += NUM_INF;
          } else {
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
            enq = qi*q_[jtype]/sqrt(r2);
            peQRealone_ += enq;

            // store new neighbor list
  //          neighCutOne_[impart].push_back(jpart);
            if ( (r2 < neighAboveSq_) && (r2 > neighBelowSq_) ) {
              neighOne_[impart].push_back(jpart);
            }
          }
        }
      }
    }
    // standard long range correction contribution of ipart
    peLRCone_ += lrcOne(ipart);
  }

  return peLJone_ + peLRCone_ + peQRealone_;
}

/**
 * Returns the total potential energy of the system
 */
double PairLJCoul::peTot() {
  return peLJ_ + peLRC_ + peQReal_;
}

/**
 * compute standard long range contributions of all particles
 */
void PairLJCoul::lrcConf() {
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
double PairLJCoul::lrcOne(
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
 *  compute interaction contributions of multiple particles
 *  flag==0, energy of old configuration (no moves, ins or dels)
 *  flag==1, new configuration, same number of particles
 *  flag==2, old configuration, preparing to delete (same as 0)
 *  flag==3, just inserted particle
 */
double PairLJCoul::multiPartEner(const vector<int> multiPart,
  const int flag) {
  multiPartEnerReal(multiPart, flag);

  double detot = 0;
  if (flag == 0 || flag == 1) {
    if (cheapEnergy_) {
      detot = peLJone_;
    } else {
      detot = peLJone_ + peLRCone_ + peQRealone_;
    }
  } else if (flag == 2) {
    if (cheapEnergy_) {
      detot = peLJone_;
    } else {
      detot = peLJone_ + peLRCone_ + peQRealone_;
    }
  } else if (flag == 3) {
    if (cheapEnergy_) {
      detot = peLJone_;
    } else {
      detot = peLJone_ + peLRCone_ + peQRealone_;
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
void PairLJCoul::update(
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
    }
  }

  if (uptypestr.compare("update") == 0) {
    if (flag == 0) {
      peLJ_ += peLJone_ - deLJ_;
      peLRC_ += peLRCone_ - deLRC_;
      peQReal_ += peQRealone_ - deQReal_;
    }
    if (flag == 2) {
      peLJ_ -= deLJ_;
      peLRC_ -= deLRC_;
      peQReal_ -= deQReal_;
    }
    if (flag == 3) {
      peLJ_ += deLJ_;
      peLRC_ += deLRC_;
      peQReal_ += deQReal_;
    }
  }
}

/**
 * initialize bulk simulation of SPCE waters, atoms listed as oxygen then two
 * hydrogens
 */
void PairLJCoul::initBulkSPCE() {
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
void PairLJCoul::initBulkSPCE(
  const double alphatmp,    //!< Ewald alpha parameter
  const int kmax) {         //!< Ewald kmax
  initBulkSPCE();
}

/**
 * check size of class variables
 */
void PairLJCoul::sizeCheck() {
}

/***
 * initialize with LAMMPS data file
 */
void PairLJCoul::initLMPData(const string fileName) {
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

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_



