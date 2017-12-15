/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./arguments.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

void Arguments::initArgs(const std::string className, const argtype &args) {
  className_ = className;
  args_ = args;
}

Arguments& Arguments::key(const std::string &key) {
  key_ = key;
  usedKeys_.push_back(key_);
  removeArg_ = false;
  return *this;
}

Arguments& Arguments::dflt(const std::string &defaultVal) {
  defaultVal_ = defaultVal;
  return *this;
}

bool Arguments::empty() {
  ASSERT(!key_.empty(), "key must be set before");
  auto pair = args_.find(key_);
  if (pair != args_.end()) {
    return false;
  }
  return true;
}

std::string Arguments::str() {
  ASSERT(!key_.empty(), "key must be set before");
  auto pair = args_.find(key_);
  std::string str;
  if (pair != args_.end()) {
    str = pair->second;

    // remove parsed arg from args_
    if (removeArg_) {
      args_.erase(pair);
      removeArg_ = false;
    }
  } else {
    ASSERT(!defaultVal_.empty(), "key(" << key_ << ") is required for args "
      << "when no default value is set.");
    str = defaultVal_;
  }
  key_ = "";
  defaultVal_ = "";
  return str;
}

bool Arguments::checkAllArgsUsed(){
  usedKeys_.push_back("");
  for (const auto &pair : args_) {
    ASSERT(findInList(pair.first, usedKeys_),
      "Key arg(" << pair.first <<") given to class(" << className_
      << ") is not recognized. All keywords provided in args must be "
      << "used, otherwise, a simple typo in keyword arguments would go "
      << "unnoticed. If you are a developer, check that you have processed "
      << "all possible keywords in the appropriate constructor. If you are "
      << "a user, check for any typos in your provided keyword arguments and "
      << "check documentation for the class(" << className_ << ")");
  }
  return true;
}

Arguments& Arguments::rm() {
  removeArg_ = true;
  return *this;
}

int Arguments::integer() {
  const std::string strVal = str();
  std::stringstream errmsg;
  int intVal = -1;
  errmsg << "Argument({" << usedKeys_.back() << ", " << strVal << "}) was "
    << "expected to be an integer";
  try {
    intVal = stoi(strVal);
  } catch (...) {
    ASSERT(0, errmsg.str());
  }
  const double dble = stod(strVal);
  ASSERT(fabs(dble - intVal) < 10*DTOL, errmsg.str());
  return intVal;
}

double Arguments::dble() {
  const std::string strVal = str();
  double dbleVal = -1;
  try {
    dbleVal = stod(strVal);
  } catch (...) {
    ASSERT(0, "Argument({" << usedKeys_.back() << ", " << strVal << "}) was "
    << "expected to be a double precision floating point number.");
  }
  return dbleVal;
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_
