
#ifndef FEASST_CORE_SITE_TEST_H_
#define FEASST_CORE_SITE_TEST_H_

#include "system/test/position_test.h"
#include "configuration/include/site.h"

namespace feasst {

inline Site default_site() {
  Site site;
  site.set_position(default_position());
  site.set_type(0);
  return site;
}

}  // namespace feasst

#endif  // FEASST_CORE_SITE_TEST_H_
