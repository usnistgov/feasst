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
#include "atom.h"

TEST(Atom, mol) {
  auto space = feasst::makeSpace({{"dimen", "3"}});
  space->addMolInit("../forcefield/data.spce");
  for (int i = 0; i < 10; ++i) space->addMol(0);
  space->addMolInit("../forcefield/data.tip3p");
  for (int i = 0; i < 10; ++i) space->addMol(1);
  feasst::Atom atom;
  space->addMol(1);
  atom.addPart(*(space.get()));
  EXPECT_EQ(21*3, atom.size());
  EXPECT_EQ(0, atom.intVal(50));
  EXPECT_EQ(0, atom.intVal(10*3+0));
  EXPECT_EQ(0, atom.intVal(10*3-1));

  // delete particle imol 1
  vector<int> mpart = space->imol2mpart(1);
  space->delPart(mpart);
  for (int iPart = 0; iPart < static_cast<int>(mpart.size()); ++iPart) {
    atom.delPart(iPart);
  }
  EXPECT_EQ(20*3, atom.size());
  EXPECT_EQ(0, atom.intVal(50));
  EXPECT_EQ(0, atom.intVal(9*3+0));
  EXPECT_EQ(0, atom.intVal(9*3-1));

  // add one of each
  space->addMol(0);
  space->addMol(1);
  atom.addPart(*(space.get()));
  EXPECT_EQ(22*3, atom.size());
  EXPECT_EQ(0, atom.intVal(21*3-1));
  EXPECT_EQ(0, atom.intVal(22*3-1));

  // delete imol 15 with "fastDel" style (e.g., swap before del)
  mpart = space->imol2mpart(15);
  for (int i = static_cast<int>(mpart.size()) - 1; i >= 0; --i) {
    const int ipart = mpart[i];
    atom.swap(ipart, space->natom() - 1);
  }
  for (int i = static_cast<int>(mpart.size()) - 1; i >= 0; --i) {
    atom.delPart(space->natom() - mpart.size() + i);
  }
  space->delPart(mpart);
  EXPECT_EQ(21*3, atom.size());
  EXPECT_EQ(0, atom.intVal(50));
  EXPECT_EQ(0, atom.intVal(9*3+0));
  EXPECT_EQ(0, atom.intVal(9*3-1));
  EXPECT_EQ(0, atom.intVal(21*3-1));

  // test restarting
  atom.writeRestart("tmp/atomrst");
  feasst::Atom atom2("tmp/atomrst");
  EXPECT_EQ(21*3, atom2.size());
  EXPECT_EQ(0, atom2.intVal(50));
  EXPECT_EQ(0, atom2.intVal(9*3+0));
  EXPECT_EQ(0, atom2.intVal(9*3-1));
  EXPECT_EQ(0, atom2.intVal(21*3-1));
}

