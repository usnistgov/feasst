#include <gtest/gtest.h>
#include "pair_wall.h"

TEST(PairWall, wall) {
  feasst::Space space(3,0);
  space.lset(8, 0);
  space.lset(30, 1);
  space.lset(30, 2);
  std::string addMolType("../forcefield/data.lj");
  space.addMolInit(addMolType);
  for (int iMol = 0; iMol < 50; ++iMol) {
    space.addMol(addMolType);
  }
  space.lset(30, 0);
  
  feasst::Barrier barrier;
  barrier.addOrthogonalPlanar(5, 1, 0);
  barrier.addOrthogonalPlanar(-5, -1, 0);
  
  feasst::PairWall pair(&space, &barrier);
  pair.initEnergy();
  EXPECT_NEAR(0, pair.peTot(), 1e-14);
  pair.printxyz("tmp/wall5", 1);
  
  for (int iMol = 0; iMol < 50; ++iMol) {
    space.addMol(addMolType);
  }
  pair.initEnergy();
  EXPECT_LT(std::numeric_limits<double>::max()/1e10, pair.peTot());
}

