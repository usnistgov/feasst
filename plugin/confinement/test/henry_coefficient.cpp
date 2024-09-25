#include "utils/test/utils.h"
#include "confinement/include/henry_coefficient.h"

namespace feasst {

TEST(HenryCoefficient, serialize) {
  HenryCoefficient obj;
  HenryCoefficient obj2 = test_serialize(obj);
}

}  // namespace feasst
