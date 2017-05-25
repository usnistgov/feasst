#include "random.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

Random::Random(const unsigned long long iseed
  ) : seed_(iseed) {
  verbose_ = 0;
  className_.assign("Random");
}

Random::Random(const char* fileName) {
  verbose_ = 0;
  className_.assign("Random");
  (void) fileName;  // avoid unused parameter warning
}

void Random::seed(const unsigned long long iseed) {
  seed_ = iseed;
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_



