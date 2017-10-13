/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#include "./pair_ss.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

PairSS::PairSS(Space* space, const double rCut)
  : Pair(space, rCut) {
  defaultConstruction_();
}

PairSS::PairSS(Space* space, const char* fileName)
  : Pair(space, fileName) {
  defaultConstruction_();
  string str = fstos("npotentialparameter", fileName);
  if (!str.empty()) {
    initExponent(stoi(str));
  }
}

void PairSS::defaultConstruction_() {
  className_.assign("PairSS");
  initExponent();
}

void PairSS::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  file << "npotentialparameter " << n_ << endl;
}

void PairSS::initEnergy() {
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
          sigij = sigij_[type[ipart]][type[jpart]];
          peTot_ += pow(sigij*sigij/r2, n_/2.);
// std::streamsize ss = cout.precision();
// cout << std::setprecision(std::numeric_limits<long double>::digits10+2)
//   << "peSR " << peTot_ << " for ipart " << ipart << " jpart " << jpart
//   << " eps " << epsij_[type[ipart]][type[jpart]] << " sig "
//   << sigij_[type[ipart]][type[jpart]] << " r2 " << r2 << endl;
// cout << std::setprecision(ss);
        }
      }
    }
  }
}

double PairSS::multiPartEner(
  const vector<int> mpart,
  const int flag) {
  if (flag == 0) {}  // remove unused parameter warning

  // zero potential energy contribution of particle ipart
  peSRone_ = 0;
  return multiPartEnerAtomCut(mpart);
}

void PairSS::multiPartEnerAtomCutInner(const double &r2,
  const int &itype,
  const int &jtype) {
  const double sigij = sigij_[itype][jtype];
  peTot_ += pow(sigij*sigij/r2, n_/2.);
  // cout << "in Inner: r2 " << r2 << " itype " << itype << " jtype " << jtype
  //   << " peSRone " << peSRone_ << " sigij " << sigij << endl;
}

void PairSS::update(const vector<int> mpart,    //!< particles involved in move
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

double PairSS::b2reduced() {
  ASSERT(n_ == 12, "virial coefficient hard coded for n=12 for gamma[1-3/n]");
  const double gamma = 1.225416702465178;
  const double sigij = sigij_[0][0];
  return 2.*PI/3.*pow(sigij, 3)*gamma;
}

shared_ptr<PairSS> makePairSS(Space* space, const double rCut) {
  return make_shared<PairSS>(space, rCut);
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_




