/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./pair_tabular_1d.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

PairTabular1D::PairTabular1D(Space* space)
  : Pair(space, 0.) {
  defaultConstruction_();
}

PairTabular1D::PairTabular1D(Space* space,
  const char* fileName)
  : Pair(space, fileName) {
  defaultConstruction_();
  string tabfile = fstos("tabFileName", fileName);
  readTable(tabfile.c_str());
}

void PairTabular1D::defaultConstruction_() {
  className_.assign("PairTabular1D");
  tol_ = 0.75;
}

void PairTabular1D::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  file << "# tabFileName " << tabFileName_ << endl;
}

void PairTabular1D::pairSiteSite_(const int &iSiteType, const int &jSiteType,
  double * energy, double * force, int * neighbor, const double &dx,
  const double &dy, const double &dz) {
  const double r2 = dx*dx + dy*dy + dz*dz;
  const double rCutInner = rCutInner_[iSiteType][jSiteType];
  if (r2 < rCutInner*rCutInner) {
    *energy = NUM_INF;
  } else {
    *energy = peTable_[iSiteType][jSiteType]->interpolate(r2);
  }
  *neighbor = 1;
}

void PairTabular1D::readTable(const char* fileName) {
  tabFileName_.assign(fileName);

  // resize rCutInner and tabij
  rCutInner_.resize(space_->nParticleTypes(), vector<double>(
    space_->nParticleTypes()));
  peTable_.resize(space_->nParticleTypes(), vector<shared_ptr<Table> >(
    space_->nParticleTypes()));

  // loop through each unique pair of particle types
  for (int iType = 0; iType < space_->nParticleTypes(); ++iType) {
    for (int jType = iType; jType < space_->nParticleTypes(); ++jType) {
      // initialize table file name
      stringstream ss;
      ss << fileName << "i" << iType << "j" << jType;

      // read rCutInner and rCut
      rCutInner_[iType][jType] = sqrt(fstod("dimmin0", ss.str().c_str()));
      rCutInner_[jType][iType] = rCutInner_[iType][jType];
      rCutijset(iType, jType, sqrt(fstod("dimmax0", ss.str().c_str())));

      // make the tables
      peTable_[iType][jType] = make_shared<Table>(ss.str().c_str());
      peTable_[jType][iType] = peTable_[iType][jType];
    }
  }
}

void PairTabular1D::setInterpolator(const char* name) {
  for (int itype = 0; itype < static_cast<int>(peTable_.size()); ++itype) {
    for (int jtype = 0; jtype < static_cast<int>(peTable_[itype].size());
         ++jtype) {
      peTable_[itype][jtype]->setInterpolator(name);
    }
  }
}

PairTabular1D* PairTabular1D::clone(Space* space) const {
  PairTabular1D* p = new PairTabular1D(*this);
  p->reconstruct(space);
  return p;
}

shared_ptr<PairTabular1D> makePairTabular1D(Space* space) {
  return make_shared<PairTabular1D>(space);
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

