/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./pair_field.h"

namespace feasst {

PairField::PairField(Space* space) : Pair(space) {
  defaultConstruction_();
}

PairField::PairField(Space* space, const char* fileName) : Pair(space, fileName) {
  defaultConstruction_();

  deWall_ = fstod("deWall", fileName);
}

void PairField::defaultConstruction_() {
  className_.assign("PairField");
}

void PairField::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  file << "# deWall " << deWall_ << endl;
}

double PairField::multiPartEner(const vector<int> mpart, const int flag) {
  if (flag == 0) {}; //remove unused parameter warning
  return fieldLoopSite_(mpart);
}

double PairField::fieldLoopSite_(const vector<int> &siteList) {
  const vector<double> &x = space_->x();
  peSRone_ = 0;
  double energy, force;
  double xi, yi = 0., zi = 0.;
  for (int impart = 0; impart < int(siteList.size()); ++impart) {
    const int iPart = siteList[impart];
    const int iType = space_->type()[iPart];
    xi = x[iPart*dimen_];
    if (dimen_ > 1) {
      yi = x[iPart*dimen_ + 1];
    }
    if (dimen_ > 2) {
      zi = x[iPart*dimen_ + 2];
    }
    fieldSite_(iType, &energy, &force, xi, yi, zi);
    peSRone_ += energy;
  }
  return peSRone_;
}

void PairField::initEnergy() {
  std::fill(pe_.begin(), pe_.end(), 0.);
  fill(0., f_);
  fill(0., vr_);
  // compute energy
  multiPartEner(space_->listAtoms());
  peTot_  = peSRone_;
}

void PairField::update(const vector<int> mpart,    //!< particles involved in move
                    const int flag,         //!< type of move
                    const char* uptype    //!< description of update type
  ) {
  std::string uptypestr(uptype);
  if (uptypestr.compare("store") == 0) {
    if (flag == 0 || flag == 2 || flag == 3) {
      deWall_ = peSRone_;
    }
  }

  if (uptypestr.compare("update") == 0) {
    if (flag == 0) {
      peTot_ += peSRone_ - deWall_;
    }
    if (flag == 2) {
      peTot_ -= deWall_;
    }
    if (flag == 3) {
      peTot_ += deWall_;
    }
  }
}

shared_ptr<PairField> makePairField(std::shared_ptr<Space> space) {
  return make_shared<PairField>(space.get());
}

}  // namespace feasst
