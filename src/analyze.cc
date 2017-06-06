#include "./analyze.h"

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
}

/**
 * defaults in constructor
 */
void Analyze::defaultConstruction() {
  nFreq_ = 1;
  nFreqPrint_ = 1;
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
 * write restart file
 */
void Analyze::writeRestartBase(const char* fileName) {
  fileBackUp(fileName);
  std::ofstream file(fileName);
  file << "# className " << className_ << endl;
  file << "# nFrequency " << nFreq_ << endl;
  file << "# nPrintFrequency " << nFreqPrint_ << endl;
  if (!fileName_.empty()) file << "# fileName " << fileName_ << endl;
}

/**
 *
 */
Analyze* Analyze::clone(Space* space, Pair* pair) const {
  ASSERT(space == NULL && pair == NULL, "empty method");
  return NULL;
}
shared_ptr<Analyze> Analyze::cloneImpl(Space* space, Pair* pair) const {
  ASSERT(space == NULL && pair == NULL, "empty method");
  shared_ptr<Analyze> empty;
  return empty;
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

