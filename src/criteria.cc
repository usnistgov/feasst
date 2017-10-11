/**
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of this agreement (see LICENSE.txt) and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#include "./criteria.h"
#include "./space.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

Criteria::Criteria(const double beta, const double activ)
  : beta_(beta),
    activ_(activ) {
  defaultConstruction_();
}
Criteria::Criteria(const char* fileName) {
  defaultConstruction_();
  ASSERT(fileExists(fileName),
    "restart file(" << fileName << ") doesn't exist");
  beta_ = fstod("beta", fileName);
  activ_ = exp(fstod("lnz", fileName));
  string strtmp = fstos("pressure", fileName);
  if (!strtmp.empty()) {
    pressureFlag_ = 1;
    pressure_ = stod(strtmp);
  }

  strtmp = fstos("nActivs", fileName);
  if (strtmp.empty()) {
    activVec_[0] = activ_;
  } else {
    activVec_.resize(stoi(strtmp));
    for (int ia = 0; ia < static_cast<int>(activVec_.size()); ++ia) {
      stringstream ss;
      ss << "activ" << ia;
      activVec_[ia] = fstod(ss.str().c_str(), fileName);
    }
  }

  // initialize random number generator
  initRNG(fileName);
}

void Criteria::defaultConstruction_() {
  verbose_ = 0;
  printBeta_ = 0;
  printPressure_ = 0;
  pressureFlag_ = 0;
  className_.assign("Criteria");
  activVec_.push_back(activ_);
}

void Criteria::writeRestartBase(const char* fileName) {
  fileBackUp(fileName);
  std::ofstream file(fileName);
  file << "# className " << className_ << endl
       << "# lnz " << std::setprecision(10) << log(activ_) << endl
       << "# beta " << std::setprecision(10) << beta_ << endl;
  if (pressureFlag_ == 1) file << "# pressure " << pressure_ << endl;

  if (activVec_.size() > 1) {
    file << "# nActivs " << activVec_.size() << endl;
    for (int ia = 0; ia < static_cast<int>(activVec_.size()); ++ia) {
      file << "# activ" << ia << " " << std::setprecision(10) << activVec_[ia] << endl;
    }
  }

  // write random number generator state
  writeRngRestart(fileName);
}

void Criteria::store(Pair* pair) {
  if (pair == NULL) {}
}

double Criteria::activ() const {
  ASSERT(activVec_.size() == 1, "accessing activity is ambiguous when "
    << activVec_.size() <<" activities are present");
  return activVec_[0];
}

double Criteria::activ(const int type) const {
  ASSERT(type < static_cast<int>(activVec_.size()), "Tried to access an "
    << "activity for particle type(" << type << "), but activVec size is "
    << "only " << activVec_.size() << ", therefore you likely need to use "
    << "criteria.addActiv( ) to input the activies of all particles.");
  return activVec_[type];
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_
