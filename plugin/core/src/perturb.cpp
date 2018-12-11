
#include "core/include/perturb.h"

namespace feasst {

void Perturb::store_old(System * system) {
  system_ = system;
  if (optimization_ == 0) {
    system_old_ = *system;
  }
}



}  // namespace feasst
