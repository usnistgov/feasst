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
#include "group.h"

TEST(Group, mol) {
  auto space = feasst::makeSpace({{"dimen", "3"}});
  space->addMolInit("../forcefield/data.spce");
  for (int i = 0; i < 10; ++i) space->addMol(0);
  space->addMolInit("../forcefield/data.tip3p");
  for (int i = 0; i < 10; ++i) space->addMol(1);
  feasst::Group g1;
  auto g2 = make_shared<feasst::Group>();
  g1.molid(1);
  g2->molid(1);
  space->initGroup(g2);
  space->addMol(1);
  g1.addPart(*space.get());
  EXPECT_EQ(21*3, g1.group().size());
  EXPECT_EQ(1, g1.group(10*3+0));
  EXPECT_EQ(0, g1.group(10*3-1));
  EXPECT_EQ(1, g1.group(50));
  EXPECT_EQ(1, g1.group(62));
  EXPECT_EQ(21*3, g2->group().size());
  EXPECT_EQ(1, g2->group(10*3+0));
  EXPECT_EQ(0, g2->group(10*3-1));
  EXPECT_EQ(1, g2->group(50));
  EXPECT_EQ(1, g2->group(62));

  // delete particle imol 1
  vector<int> mpart = space->imol2mpart(1);
  space->delPart(mpart);
  for (int iPart = 0; iPart < static_cast<int>(mpart.size()); ++iPart) {
    g1.delPart(iPart);
  }
  EXPECT_EQ(20*3, g1.group().size());
  EXPECT_EQ(1, g1.group(9*3+0));
  EXPECT_EQ(0, g1.group(9*3-1));
  EXPECT_EQ(1, g1.group(50));
  EXPECT_EQ(1, g1.group(59));
  EXPECT_EQ(20*3, g2->group().size());
  EXPECT_EQ(1, g2->group(9*3+0));
  EXPECT_EQ(0, g2->group(9*3-1));
  EXPECT_EQ(1, g2->group(50));
  EXPECT_EQ(1, g2->group(59));

  // add one of each
  space->addMol(0);
  space->addMol(1);
  g1.addPart(*space.get());
  EXPECT_EQ(22*3, g1.group().size());
  EXPECT_EQ(0, g1.group(21*3-1));
  EXPECT_EQ(1, g1.group(22*3-1));
  EXPECT_EQ(22*3, g2->group().size());
  EXPECT_EQ(0, g2->group(21*3-1));
  EXPECT_EQ(1, g2->group(22*3-1));

  // delete imol 15 with "fastDel" style (e.g., swap before del)
  mpart = space->imol2mpart(15);
  for (int i = static_cast<int>(mpart.size()) - 1; i >= 0; --i) {
    const int ipart = mpart[i];
    g1.swap(ipart, space->natom() - 1);
  }
  for (int i = static_cast<int>(mpart.size()) - 1; i >= 0; --i) {
    g1.delPart(space->natom() - mpart.size() + i);
  }
  space->delPart(mpart);

  EXPECT_EQ(21*3, g1.group().size());
  EXPECT_EQ(1, g1.group(50));
  EXPECT_EQ(1, g1.group(9*3+0));
  EXPECT_EQ(0, g1.group(9*3-1));
  EXPECT_EQ(0, g1.group(21*3-1));
  EXPECT_EQ(21*3, g2->group().size());
  EXPECT_EQ(1, g2->group(50));
  EXPECT_EQ(1, g2->group(9*3+0));
  EXPECT_EQ(0, g2->group(9*3-1));
  EXPECT_EQ(0, g2->group(21*3-1));

//  cout << "g " << feasst::vec2str(g1.intVal()) << endl;

  // write/read restart
  g1.writeRestart("tmp/grprst");
  feasst::Group g11("tmp/grprst");
  EXPECT_EQ(21*3, g11.group().size());
  EXPECT_EQ(1, g11.group(50));
  EXPECT_EQ(1, g11.group(9*3+0));
  EXPECT_EQ(0, g11.group(9*3-1));
  EXPECT_EQ(0, g11.group(21*3-1));

  space->writeRestart("tmp/grprst2");
  feasst::Space space2("tmp/grprst2");
  shared_ptr<feasst::Group> g22 = space2.groups()[0];
  EXPECT_EQ(21*3, g22->group().size());
  EXPECT_EQ(1, g22->group(50));
  EXPECT_EQ(1, g22->group(9*3+0));
  EXPECT_EQ(0, g22->group(9*3-1));
  EXPECT_EQ(0, g22->group(21*3-1));

}

