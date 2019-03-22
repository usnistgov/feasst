
#include "flat_histogram/include/macrostate.h"

namespace feasst {

bool Macrostate::is_in_range(const System* system, const Criteria* criteria) {
  const double val = value(system, criteria);
  if (val <= histogram_.max() &&
      val >= histogram_.min()) {
    return true;
  }
  return false;
}

}  // namespace feasst
