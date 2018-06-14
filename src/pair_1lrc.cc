/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./pair_1lrc.h"

namespace feasst {

PairLRC::PairLRC(Space* space, const argtype &args)
  : Pair(space, args) {
  defaultConstruction_();
}

PairLRC::PairLRC(Space* space, const char* fileName)
  : Pair(space, fileName) {
  defaultConstruction_();
}

void PairLRC::defaultConstruction_() {
  className_.assign("PairLRC");
}

void PairLRC::initLRC() {
  lrcFlag = 1;
  lrcPreCalc_.resize(epsij_.size(), vector<double>(epsij_.size()));
  ASSERT(sigrefFlag_ == 0, "sigref not implemented here");
  for (unsigned int i = 0; i < epsij_.size(); ++i) {
    for (unsigned int j = 0; j < epsij_.size(); ++j) {
      const double sig = sigij_[i][j];
      const double eps = epsij_[i][j];
      double rc = 0.;
      if (rCutij_.size() == 0 || atomCut_ == 0) {
        rc = rCut_;
      } else {
        rc = rCutij_[i][j];
        // cout << "rc" << i << j << " " << rc << endl;
      }
      if ( (fabs(sig) > DTOL) && (fabs(eps) > DTOL) ) {
        lrcPreCalc_[i][j] = (8./3.)*PI*eps*pow(sig, 3)*(pow(sig/rc, 9)/3.
          - pow(sig/rc, 3));
        lrcPreCalc_[j][i] = lrcPreCalc_[i][j];
      }
    }
  }
//  cout << "Init LRCs" << endl << vec2str(lrcPreCalc_) << endl;
}

double PairLRC::computeLRC(const vector<int> msite) {
  double enlrc = 0.;
  if (lrcFlag == 0) return enlrc;
  if (lrcPreCalc_.size() == 0) {
    initLRC();
  }
  // If msite is empty, loop over all sites
  if (msite.size() == 0) {
    for (int iType = 0; iType < space_->nParticleTypes(); ++iType) {
      for (int jType = 0; jType < space_->nParticleTypes(); ++jType) {
        const int nanb = space_->nType()[iType]*space_->nType()[jType];
        enlrc += static_cast<double>(nanb)/space_->volume()
                *lrcPreCalc_[iType][jType];
      }
    }
  } else {
    // Find the number of particles of each type
    vector<int> msiteTypes(space_->nParticleTypes(), 0);
    for (int imsite = 0; imsite < static_cast<int>(msite.size()); ++imsite) {
      ++msiteTypes[space()->type()[msite[imsite]]];
    }

    for (int iType = 0; iType < space_->nParticleTypes(); ++iType) {
      for (int jType = 0; jType < space_->nParticleTypes(); ++jType) {
        const double ni = space_->nType()[iType],
                     nia = msiteTypes[iType],
                     nj = space_->nType()[jType],
                     nja = msiteTypes[jType];
        enlrc += (ni*nja + nia*nj - nia*nja)/space_->volume()
                 *lrcPreCalc_[iType][jType];
      }
    }
  }
  return enlrc;
}

}  // namespace feasst
