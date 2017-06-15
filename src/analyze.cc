#include "./analyze.h"
#include "./analyze_extensive_moments.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

Analyze::Analyze(
  Space *space,
  Pair *pair)
  : space_(space),
    pair_(pair) {
  className_.assign("Analyze");
  defaultConstruction();
}

Analyze::Analyze(Space *space,  Pair *pair, const char* fileName)
  : space_(space),
    pair_(pair) {
  className_.assign("Analyze");

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

/**
 * defaults in constructor
 */
void Analyze::defaultConstruction() {
  nFreq_ = 1;
  nFreqPrint_ = 1;
  production_ = 1;
}

/**
 * reset object pointers
 */
void Analyze::reconstruct(Space* space, Pair *pair) {
  space_ = space;
  pair_ = pair;
  Base::reconstruct();
}

/**
 * factory method
 */
shared_ptr<Analyze> Analyze::makeAnalyze(
  Space* space,
  Pair* pair,
  const char* fileName) {
  ASSERT(fileExists(fileName),
         "restart file(" << fileName << ") doesn't exist");
  string anatypestr = fstos("className", fileName);
  shared_ptr<Analyze> ana;
  if (anatypestr.compare("AnalyzeExtensiveMoments") == 0) {
    ana = make_shared<AnalyzeExtensiveMoments>(space, pair, fileName);
  } else {
    ASSERT(0, "unrecognized analyze class(" << anatypestr << ")");
  }
  return ana;
}

/**
 * write restart file
 */
void Analyze::writeRestartBase(const char* fileName) {
  fileBackUp(fileName);
  std::ofstream file(fileName);
  file << "# className " << className_ << endl;
  file << "# nFrequency " << nFreq_ << endl;
  file << "# nPrintFrequency " << nFreqPrint_ << endl;
  file << "# production " << production_ << endl;
  if (!fileName_.empty()) file << "# fileName " << fileName_ << endl;
}

/**
 *
 */
Analyze* Analyze::clone(Space* space, Pair* pair) const {
  ASSERT(0, "base class not implemented correctly");
  return NULL;
}

shared_ptr<Analyze> Analyze::cloneImpl(Space* space, Pair* pair) const {
  ASSERT(0, "base class not implemented correctly");
  shared_ptr<Analyze> an = make_shared<Analyze>(*this);
  an->reconstruct(space, pair);
  return an;
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_
