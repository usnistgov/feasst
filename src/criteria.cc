#include "./criteria.h"
#include "./criteria_metropolis.h"
#include "./criteria_wltmmc.h"
#include "./criteria_mayer.h"
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
       << "# lnz " << log(activ_) << endl
       << "# beta " << beta_ << endl;
  if (pressureFlag_ == 1) file << "# pressure " << pressure_ << endl;

  if (activVec_.size() > 1) {
    file << "# nActivs " << activVec_.size() << endl;
    for (int ia = 0; ia < static_cast<int>(activVec_.size()); ++ia) {
      file << "# activ" << ia << " " << activVec_[ia] << endl;
    }
  }

  // write random number generator state
  writeRngRestart(fileName);
}

void Criteria::store(const Space* space, Pair* pair) {
  if (space == NULL) if (pair == NULL) {}
}

Criteria* Criteria::makeCriteria(const char* fileName
  ) {
  ASSERT(fileExists(fileName),
    "restart file(" << fileName << ") doesn't exist");
  string typestr = fstos("className", fileName);
  Criteria* criteria = NULL;
  if ( (typestr.compare("Criteria") == 0) ||
       (typestr.compare("CriteriaMetropolis") == 0) ) {
    criteria = new CriteriaMetropolis(fileName);
  } else if (typestr.compare("CriteriaWLTMMC") == 0) {
    criteria = new CriteriaWLTMMC(fileName);
  } else if (typestr.compare("CriteriaMayer") == 0) {
    criteria = new CriteriaMayer(fileName);
  } else {
    ASSERT(0, "unrecognized criteria(" << typestr << ") in factory");
  }
  return criteria;
}
  
double Criteria::activ() const {
  ASSERT(activVec_.size() == 1, "accessing activity is ambiguous when "
    << activVec_.size() <<" activities are present");
  return activVec_[0];
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

