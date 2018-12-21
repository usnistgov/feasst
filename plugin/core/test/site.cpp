#include <gtest/gtest.h>
#include <fstream>
#include "core/test/site_test.h"

namespace feasst {

TEST(Site, getset) {
  Site site = default_site();
  EXPECT_EQ(site.position().coord(0), 0.);
  Position position;
  position.set_vector({3.5, 796.4, -45.4});
  site.set_position(position);
  EXPECT_EQ(site.position().coord(0), 3.5);
  position.set_coord(1, -857.3);
  site.set_position(position);
  EXPECT_EQ(site.position().coord(1), -857.3);
}

}  // namespace feasst
