/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include <gtest/gtest.h>
#include "space.h"
#include "bond.h"

TEST(Bond, bond) {
  auto space = feasst::makeSpace({{"dimen", "3"}});
  space->addMolInit("../forcefield/data.cg1_2patch_linear");
  auto bond = make_shared<feasst::Bond>();
  space->initAtom(bond);
  space->addMol();
  space->addMol();
  space->addMolInit("../forcefield/data.lj");
  space->addMol("../forcefield/data.lj");
  space->addMol("../forcefield/data.cg1_2patch_linear");
  vector<int> intValExp = {1, 0, 0,         // first mol
                                 1, 0, 0,   // second mol
                                 bond->defaultInt(), // lj
                                 1, 0, 0};  // last
  EXPECT_EQ(bond->intVal(), intValExp);
  const vector<int> mpart = space->imol2mpart(1);
  space->delPart(mpart);
  intValExp = {1, 0, 0,      // first mol
               1, 0, 0,      // second mol
               bond->defaultInt()};   // lj
  EXPECT_EQ(bond->intVal(), intValExp);
}
