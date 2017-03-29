/**
 * \file
 *
 * \brief trial moves for Monte Carlo
 */

#include "./criteria.h"
#include "./criteria_metropolis.h"
#include "./criteria_wltmmc.h"
#include "./space.h"

/**
 * Constructor
 */
Criteria::Criteria(const double beta,  //!< inverse temperature
  const double activ)  //!< activity
  : beta_(beta),
    activ_(activ) {
  defaultConstruction();
}
Criteria::Criteria(const char* fileName) {
  defaultConstruction();
  ASSERT(myFileExists(fileName),
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
    activVec[0] = activ_;
  } else {
    activVec.resize(stoi(strtmp));
    for (int ia = 0; ia < static_cast<int>(activVec.size()); ++ia) {
      stringstream ss;
      ss << "activ" << ia;
      activVec[ia] = fstod(ss.str().c_str(), fileName);
    }
  }
  
  // initialize random number generator
  initRNG(fileName);
}

/**
 * defaults in constructor
 */
void Criteria::defaultConstruction() {
  verbose_ = 0;
  printBeta_ = 0;
  printPressure_ = 0;
  pressureFlag_ = 0;
  className_.assign("Criteria");
  activVec.push_back(activ_);
}

/**
 * write restart file
 */
void Criteria::writeRestartBase(const char* fileName) {
  fileBackUp(fileName);
  std::ofstream file(fileName);
  file << "# className " << className_ << endl
       << "# lnz " << log(activ_) << endl
       << "# beta " << beta_ << endl;
  if (pressureFlag_ == 1) file << "# pressure " << pressure_ << endl;

  if (activVec.size() > 1) {
    file << "# nActivs " << activVec.size() << endl;
    for (int ia = 0; ia < static_cast<int>(activVec.size()); ++ia) {
      file << "# activ" << ia << " " << activVec[ia] << endl;
    }
  }

  // write random number generator state
  writeRngRestart(fileName);
}

/**
 * store macrostate variables of old configuration
 */
void Criteria::store(const Space* space, Pair* pair) {
  if (space == NULL) if (pair == NULL) {}
}

/**
 * factory method
 */
Criteria* Criteria::makeCriteria(const char* fileName
  ) {
  ASSERT(myFileExists(fileName),
    "restart file(" << fileName << ") doesn't exist");
  string typestr = fstos("className", fileName);
  Criteria* criteria = NULL;
  if ( (typestr.compare("Criteria") == 0) ||
       (typestr.compare("CriteriaMetropolis") == 0) ) {
    criteria = new CriteriaMetropolis(fileName);
  } else if (typestr.compare("CriteriaWLTMMC") == 0) {
    criteria = new CriteriaWLTMMC(fileName);
  } else {
    ASSERT(0, "unrecognized criteria(" << typestr << ") in factory");
  }
  return criteria;
}
  
/**
 * external access to activity
 */
double Criteria::activ() const {
  ASSERT(activVec.size() == 1, "accessing activity is ambiguous when "
    << activVec.size() <<" activities are present");
  //return activ_;
  return activVec[0];
}
