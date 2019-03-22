
#include <sstream>
#include "core/include/accumulator.h"
#include "core/include/debug.h"
#include "flat_histogram/include/bias.h"

namespace feasst {

std::string Bias::write_per_bin(const int bin) const {
  std::stringstream ss;
  ss << ln_macro_prob().value(bin);
  return ss.str();
}

int Bias::bin_(
    const int macrostate_old,
    const int macrostate_new,
    const bool is_accepted) {
  int bin = -1;
  if (is_accepted) {
    bin = macrostate_new;
  } else {
    bin = macrostate_old;
  }
  DEBUG("bin " << bin);
  return bin;
}

}  // namespace feasst
