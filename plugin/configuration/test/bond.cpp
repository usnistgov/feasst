#include "utils/test/utils.h"
#include "configuration/include/bond.h"

namespace feasst {

TEST(Bond, serialize) {
  Bond bond;
  bond.add_site_index(12);
  bond.add_site_index(1);
  bond.set_type(45);
  bond.add_property("banana", 12);

  Bond bond2 = test_serialize(bond, "847 1 1 banana 1 12 Bond 743 45 2 12 1 0 ");
  EXPECT_EQ(bond2.type(), 45);
  EXPECT_EQ(bond2.property("banana"), 12.0);
  EXPECT_EQ(bond2.site(0), 12);
}

TEST(Angle, serialize) {
  Angle angle;
  angle.add_site_index(12);
  angle.add_site_index(1);
  angle.set_type(45);
  angle.add_property("banana", 12);

  Angle angle2 = test_serialize(angle, "847 1 1 banana 1 12 Angle 743 45 2 12 1 0 ");
  EXPECT_EQ(angle2.class_name(), "Angle");
}

}  // namespace feasst
