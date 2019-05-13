#include <fstream>
#include "utils/test/utils.h"
#include "configuration/test/site_test.h"

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

  // serialize
  site.set_type(125);
  site.add_property("banana", 12);
  Site site2 = test_serialize(site);
  EXPECT_EQ(site.type(), site2.type());
  EXPECT_EQ(site.properties().str(), site2.properties().str());
  EXPECT_EQ(site.position().coord(), site2.position().coord());
  EXPECT_EQ(site.is_director(), site2.is_director());
}

}  // namespace feasst
