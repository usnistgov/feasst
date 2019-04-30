#include <gtest/gtest.h>
#include "core/include/bond.h"

namespace feasst {

TEST(Bond, serialize) {
  Bond bond;
  bond.set_name("custom_bond");
  bond.add_site_index(12);
  bond.add_site_index(1);
  bond.set_type(45);
  bond.add_property("banana", 12);

  std::stringstream ss;
  bond.serialize(ss);
  Bond bond2(ss);
  EXPECT_EQ(bond2.type(), 45);
  EXPECT_EQ(bond2.property("banana"), 12.0);
  EXPECT_EQ(bond2.site(0), 12);
  EXPECT_EQ("1 1 banana 1 12 1 45 1 custom_bond 2 12 1 ", ss.str());
}

TEST(Angle, serialize) {
  Angle angle;
  angle.add_site_index(12);
  angle.add_site_index(1);
  angle.set_type(45);
  angle.add_property("banana", 12);

  std::stringstream ss;
  angle.serialize(ss);
  Angle angle2(ss);
  EXPECT_EQ(angle2.name(), "angle");
  EXPECT_EQ("1 1 banana 1 12 1 45 1 angle 2 12 1 ", ss.str());
}

}  // namespace feasst
