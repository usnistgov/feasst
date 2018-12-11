#include <gtest/gtest.h>
#include "ewald/include/ewald.h"

TEST(Ewald, ewald) {
  feasst::Ewald ewald;
  EXPECT_FALSE(ewald.enabled());
}
