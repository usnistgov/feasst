#include "prefetch/include/pool.h"

namespace feasst {

const std::string Pool::str() const {
  std::stringstream ss;
  ss << index_ << " " << ln_prob_ << " " << accepted_;
  return ss.str();
}

}  // namespace feasst
