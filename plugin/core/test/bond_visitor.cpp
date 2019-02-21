#include <gtest/gtest.h>
#include "core/include/bond_visitor.h"
#include "core/test/system_test.h"

namespace feasst {

TEST(BondVisitor, spce) {
  Configuration config = spce_sample();
  EXPECT_EQ(2, config.particle_type(0).num_bonds());
  EXPECT_EQ(1, config.unique_type(0).num_bonds());
  EXPECT_NEAR(1., config.unique_type(0).bond(0).property("l0"), NEAR_ZERO);
  BondVisitor visitor;
  BondSquareWell model;
  visitor.compute(model, config);
  EXPECT_NEAR(15*NEAR_INFINITY, visitor.energy(), NEAR_INFINITY/1e10);
  AngleSquareWell angle;
  visitor.compute(angle, config);
  EXPECT_NEAR(14*NEAR_INFINITY, visitor.energy(), NEAR_INFINITY/1e10);
}

}  // namespace feasst
