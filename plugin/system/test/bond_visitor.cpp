#include "utils/test/utils.h"
#include "configuration/include/utils.h"
#include "system/include/bond_visitor.h"
#include "system/include/bond_square_well.h"
#include "system/include/angle_square_well.h"

namespace feasst {

TEST(BondVisitor, spce) {
  Configuration config = spce_sample1();
  EXPECT_EQ(2, config.particle_type(0).num_bonds());
  EXPECT_EQ(1, config.unique_type(0).num_bonds());
  EXPECT_NEAR(1., config.unique_type(0).bond(0).property("length"), NEAR_ZERO);
  BondVisitor visitor;
  BondSquareWell model;
  visitor.compute_all(config);
  DEBUG(visitor.energy()/NEAR_INFINITY);
  EXPECT_EQ(0, visitor.energy());
  AngleSquareWell angle;
  visitor.compute_all(config);
  EXPECT_EQ(0, visitor.energy());

  BondVisitor visitor2 = test_serialize(visitor);
  BondSquareWell model2 = test_serialize(model);
  AngleSquareWell angle2 = test_serialize(angle);
}

}  // namespace feasst
