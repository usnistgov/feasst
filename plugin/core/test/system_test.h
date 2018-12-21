
#ifndef FEASST_CORE_SYSTEM_TEST_H_
#define FEASST_CORE_SYSTEM_TEST_H_

#include <vector>
#include <utility>
#include <string>
#include "core/include/system.h"
#include "core/test/configuration_test.h"

namespace feasst {

inline System default_system() {
  System sys;
  Configuration config = default_configuration();
  sys.add_configuration(config);
  return sys;
}

}  // namespace feasst

#endif  // FEASST_CORE_SYSTEM_TEST_H_
