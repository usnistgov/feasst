
#ifndef FEASST_CORE_POSITION_TEST_H_
#define FEASST_CORE_POSITION_TEST_H_

#include "math/include/position.h"

namespace feasst {

inline Position default_position() {
  Position pos;
  pos.set_to_origin_3D();
  return pos;
}

}  // namespace feasst

#endif  // FEASST_CORE_POSITION_TEST_H_
