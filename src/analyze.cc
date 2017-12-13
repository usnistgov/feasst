/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./analyze.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

Analyze::Analyze(Pair *pair, const argtype &args)
  : pair_(pair) {
  defaultConstruction_();
  argparse_.initArgs(className_, args);

  // parse nFreq
  if (!argparse_.key("nFreq").empty()) {
    initFreq(stoi(argparse_.str()));
  }

  // parse nFreqPrint
  if (!argparse_.key("nFreqPrint").empty()) {
    initFreqPrint(stoi(argparse_.str()));
  }

}

Analyze::Analyze(Pair *pair, const char* fileName)
  : pair_(pair) {
  defaultConstruction_();

  ASSERT(fileExists(fileName),
    "restart file(" << fileName << ") doesn't exist");

  nFreq_ = fstoi("nFrequency", fileName);
  nFreqPrint_ = fstoi("nPrintFrequency", fileName);
  string strtmp = fstos("fileName", fileName);
  if (!strtmp.empty()) {
    fileName_ = strtmp;
  }
  strtmp = fstos("production", fileName);
  if (!strtmp.empty()) {
    production_ = stoi(strtmp);
  }
}

void Analyze::defaultConstruction_() {
  className_.assign("Analyze");
  initFreq();
  initPrintFreq();
  initProduction();
}

void Analyze::reconstruct(Pair *pair) {
  pair_ = pair;
  Base::reconstruct();
}

void Analyze::writeRestartBase(const char* fileName) {
  fileBackUp(fileName);
  std::ofstream file(fileName);
  file << "# className " << className_ << endl;
  file << "# nFrequency " << nFreq_ << endl;
  file << "# nPrintFrequency " << nFreqPrint_ << endl;
  file << "# production " << production_ << endl;
  if (!fileName_.empty()) file << "# fileName " << fileName_ << endl;
}

Analyze* Analyze::clone(Pair* pair) const {
  ASSERT(0, "base class not implemented correctly");
  return NULL;
}

shared_ptr<Analyze> Analyze::cloneImpl(Pair* pair) const {
  ASSERT(0, "base class not implemented correctly");
  shared_ptr<Analyze> an = make_shared<Analyze>(*this);
  an->reconstruct(pair);
  return an;
}

shared_ptr<Analyze> makeAnalyze(Pair* pair) {
  return make_shared<Analyze>(pair);
}

shared_ptr<Analyze> makeAnalyze(Pair* pair, const argtype &args) {
  return make_shared<Analyze>(pair, args);
}

shared_ptr<Analyze> makeAnalyze() {
  return make_shared<Analyze>();
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_
