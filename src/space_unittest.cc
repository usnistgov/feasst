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

using namespace feasst;

TEST(Space, readwritexyz) {
  Space s(3);
  s.init_config(12);
  EXPECT_EQ(1, s.checkSizes());
  s.printXYZ("tmp/tmp", 1);
  EXPECT_EQ(1, s.checkSizes());
  s.xset(500, 0,0);
  s.readXYZ("tmp/tmp.xyz");
  EXPECT_EQ(s.x(0,0), 0);
  EXPECT_EQ(1, s.checkSizes());
}

TEST(Space, delPartaddPart) {
  Space s(3);
  s.init_config(12);
  s.tagAtom(3); s.tagAtom(2);
  vector<int> mpart(1, 2);
  s.delPart(mpart);
  EXPECT_EQ(s.natom(), 11);
  EXPECT_EQ(s.natom(),int(s.type().size()));
  EXPECT_EQ(s.natom(),int(s.mol().size()));
  EXPECT_EQ(s.natom(),int(s.moltype().size()));
  std::stringstream ss;
  ss << s.install_dir() << "/forcefield/data.lj";
  EXPECT_EQ(0, s.moltype().front().compare(ss.str()));
  EXPECT_EQ(0, s.moltype().back().compare(ss.str()));
  EXPECT_EQ(s.natom(), s.nType()[0]);
  EXPECT_EQ(1, s.nParticleTypes());
  EXPECT_EQ(2, int(s.tag().size()));
  EXPECT_EQ(3, s.tag()[0]);
  EXPECT_EQ(-1, s.tag()[1]);
  s.tagAtomPopBack();
  EXPECT_EQ(1, int(s.tag().size()));
  EXPECT_EQ(3, s.tag()[0]);
  s.tagAtomPopBack();
  EXPECT_EQ(0, int(s.tag().size()));

  vector<double> x(s.dimen());
  x[0] = 5; x[1] = 5; x[2] = 5;
  // this type of addPart doesn't correctly keep up with moltype_ and possibly other arrays
  s.addPart(x, 1, s.mol().back()+1);
  EXPECT_EQ(s.natom(), 12);
  EXPECT_EQ(s.x(11,0), 5);
  EXPECT_EQ(s.natom(),int(s.type().size()));
  EXPECT_EQ(s.natom(),int(s.mol().size()));
  EXPECT_EQ(1,s.type().at(11));
  EXPECT_EQ(11,s.mol().at(11));
  x[0] = 6; x[1] = 6; x[2] = 6;
  s.addPart(x, 1, s.mol().back() + 1);
  EXPECT_EQ(s.natom(), 13);
  EXPECT_EQ(s.x(12,0), 6);
  EXPECT_EQ(s.natom(),int(s.type().size()));
  EXPECT_EQ(s.natom(),int(s.mol().size()));
  EXPECT_EQ(12,s.mol().at(12));
  EXPECT_EQ(1,s.type().at(12));
  EXPECT_EQ(s.natom() - 2, s.nType()[0]);
  EXPECT_EQ(2, s.nType()[1]);
  EXPECT_EQ(2, s.nParticleTypes());
  EXPECT_EQ(s.natom()-2,int(s.moltype().size()));
  EXPECT_EQ(0, s.moltype().front().compare(ss.str()));
  EXPECT_EQ(0, s.moltype().back().compare(ss.str()));
}

TEST(Space, vol) {
  int dim=3, natom=12;
  Space s(dim);
  s.init_config(natom);
  EXPECT_EQ(pow(natom*dim, 3), s.volume());
}

TEST(Space, readwater) {
  int dim=3;
  Space s(dim);
  s.readXYZBulk(3, "water", "../unittest/spce/srsw/spce_sample_config_periodic1.xyz");
  EXPECT_EQ(300, s.natom());
  EXPECT_EQ(-9.011963206010, s.x(1,1));
  EXPECT_EQ(1, s.mol().at(3));
  EXPECT_EQ(0, s.type().at(3));
  EXPECT_EQ(s.natom()/3, s.nType()[0]);
  EXPECT_EQ(s.natom()*2/3, s.nType()[1]);
  EXPECT_EQ(2, s.nParticleTypes());
  EXPECT_EQ(99, s.mol().at(s.natom()-1));
  EXPECT_EQ(s.natom()/3,int(s.moltype().size()));
  EXPECT_EQ(s.natom()/3,int(s.xMolRef().size()));
  EXPECT_EQ(s.natom()/3,int(s.qMol().size())/s.qdim());
  EXPECT_EQ(0, s.moltype().front().compare("water"));
  EXPECT_EQ(0, s.moltype().back().compare("water"));
}

TEST(Space, randDispANDrestore) {
  int dim=3;
  Space s(dim);
  s.readXYZBulk(3, "water", "../unittest/spce/srsw/spce_sample_config_periodic1.xyz");
  ranInitByDate();
  vector<int> mpart(3);
  mpart[0] = 0; mpart[1] = 1; mpart[2] = 2;
  EXPECT_EQ(-5.221309047080, s.x(0,0));
  EXPECT_EQ(-9.011963206010, s.x(1,1));
  EXPECT_EQ(-9.111281355510, s.x(2,2));
  EXPECT_EQ(-1.904156489730, s.x(3,0));
  s.xStore(mpart);
  s.xStoreMulti(mpart, -1);
  EXPECT_EQ(1, int(s.xOldMulti().size()));
  s.xStoreMulti(mpart, -2);
  s.xStoreMulti(mpart, -2);
  EXPECT_EQ(3, int(s.xOldMulti().size()));
  s.randDisp(mpart, 2);
  EXPECT_EQ(-5.221309047080, s.xold().at(0).at(0));
  EXPECT_EQ(-9.011963206010, s.xold().at(1).at(1));
  EXPECT_EQ(-9.111281355510, s.xold().at(2).at(2));
  EXPECT_NE(-5.221309047080, s.x(0,0));
  EXPECT_NE(-9.011963206010, s.x(1,1));
  EXPECT_NE(-9.111281355510, s.x(2,2));
  EXPECT_EQ(-1.904156489730, s.x(3,0));
  s.restore(mpart);
  EXPECT_EQ(-5.221309047080, s.x(0,0));
  EXPECT_EQ(-9.011963206010, s.x(1,1));
  EXPECT_EQ(-9.111281355510, s.x(2,2));
  EXPECT_EQ(-1.904156489730, s.x(3,0));
  s.xStoreMulti(mpart, 0);
  EXPECT_EQ(-5.221309047080, s.x(0,0));
  EXPECT_EQ(-9.011963206010, s.x(1,1));
  EXPECT_EQ(-9.111281355510, s.x(2,2));
  s.randDisp(mpart, 2);
  s.xStoreMulti(mpart, 1);
  EXPECT_EQ(-5.221309047080, s.x(0,0));
  EXPECT_EQ(-9.011963206010, s.x(1,1));
  EXPECT_EQ(-9.111281355510, s.x(2,2));
  s.randDisp(mpart, 2);
  s.xStoreMulti(mpart, 2);
  EXPECT_EQ(-5.221309047080, s.x(0,0));
  EXPECT_EQ(-9.011963206010, s.x(1,1));
  EXPECT_EQ(-9.111281355510, s.x(2,2));
  s.randDisp(mpart, 2);
  s.xStoreMulti(mpart, -1);
  EXPECT_EQ(1, int(s.xOldMulti().size()));
}

TEST(Space, randMol) {
  int dim=3;
  Space s(dim);
  s.readXYZBulk(3, "water", "../unittest/spce/srsw/spce_sample_config_periodic1.xyz");
  ranInitByDate();
  vector<int> mpart;
  for (int i = 0; i < 1000; ++i) {
    mpart.clear();
    mpart = s.randMol();
    EXPECT_EQ(3, int(mpart.size()));
    EXPECT_LT(mpart[0], s.natom());
    EXPECT_GT(mpart[0], -1);
    mpart.clear();
    mpart = s.randMol();
    EXPECT_EQ(3, int(mpart.size()));
    EXPECT_LT(mpart[0], s.natom());
    EXPECT_GT(mpart[0], -1);
  }
}

TEST(Space, randRotateANDRestore) {
  int dim=3;
  Space s(dim);
  s.readXYZBulk(3, "water", "../unittest/spce/srsw/spce_sample_config_periodic1.xyz");
  ranInitByDate();
  vector<int> imove(3);
  imove[0] = 0; imove[1] = 1; imove[2] = 2;
  EXPECT_EQ(-5.221309047080, s.x(0,0));
  EXPECT_EQ(-9.011963206010, s.x(1,1));
  EXPECT_EQ(-9.111281355510, s.x(2,2));
  EXPECT_EQ(-1.904156489730, s.x(3,0));
  s.xStore(imove);
  s.randRotate(imove, 2);
  EXPECT_EQ(-5.221309047080, s.xold().at(0).at(0));
  EXPECT_EQ(-9.011963206010, s.xold().at(1).at(1));
  EXPECT_EQ(-9.111281355510, s.xold().at(2).at(2));
  EXPECT_NEAR(-5.22131, s.x(0,0), 1e-5);
  EXPECT_NE(-9.011963206010, s.x(1,1));
  EXPECT_NE(-9.111281355510, s.x(2,2));
  EXPECT_EQ(-1.904156489730, s.x(3,0));
  for (int i = 1; i < 3; ++i) {
    double d = 0;
    for (int dim = 0; dim < 3; ++dim) {
      d += pow(s.x(imove[i],dim) - s.x(0,dim), 2.);
    }
    EXPECT_NEAR(1, d, 1e-11);
  }
  s.restore(imove);
  EXPECT_EQ(-5.221309047080, s.x(0,0));
  EXPECT_EQ(-9.011963206010, s.x(1,1));
  EXPECT_EQ(-9.111281355510, s.x(2,2));
  EXPECT_EQ(-1.904156489730, s.x(3,0));
}

TEST(Space, addMol) {
  const int dim = 3;
  Space s(dim);
  s.readXYZBulk(3, "water", "../unittest/spce/srsw/spce_sample_config_periodic1.xyz");
  s.addMolInit("../forcefield/data.spce");
  s.addMolInit("../forcefield/data.lj");
  ranInitByDate();
  EXPECT_EQ(0, s.moltype().front().compare("water"));
  EXPECT_EQ(300, s.natom());
  s.addMol("../forcefield/data.spce");
  EXPECT_EQ(303, s.natom());
  EXPECT_EQ(101, s.nMol());
  EXPECT_EQ(0, s.moltype().back().compare("../forcefield/data.spce"));
  //s.addMol("../forcefield/data.spce");
  s.addMol("../forcefield/data.spce");
  EXPECT_EQ(306, s.natom());
  EXPECT_EQ(102, s.nMol());
  EXPECT_EQ(0, s.moltype().back().compare("../forcefield/data.spce"));
  s.addMol("../forcefield/data.lj");
  //s.addMol("atom");
  EXPECT_EQ(307, s.natom());
  EXPECT_EQ(103, s.nMol());
  EXPECT_EQ((s.natom()-1)/3, s.nType()[0]);
  EXPECT_EQ((s.natom()-1)*2/3, s.nType()[1]);
  EXPECT_EQ(3, s.nParticleTypes());
  EXPECT_EQ(0, s.moltype().front().compare("water"));
  EXPECT_EQ(0, s.moltype().back().compare("../forcefield/data.lj"));
  EXPECT_EQ(1, s.checkSizes());
}

TEST(Space, addMolANDdelANDlastMolIDVec) {
  ranInitByDate();
  const int dim = 3;
  Space s(dim);
  s.readXYZBulk(3, "water", "../unittest/spce/srsw/spce_sample_config_periodic1.xyz");
  s.addMolInit("../forcefield/data.spce");
  vector<int> mpart = s.lastMolIDVec();
  EXPECT_EQ(297, mpart[0]);
  mpart.clear();
  EXPECT_EQ(s.nMol(), int(s.mol2part().size())-1);
  for (int imol = 0; imol < s.nMol(); ++imol) {
    EXPECT_EQ(imol*3, s.mol2part()[imol]);
  }
  EXPECT_EQ(s.natom(), s.mol2part().back());
  EXPECT_EQ(s.nMol() - 1, s.mol().back());
  for (int i=0;i<s.natom();++i) EXPECT_EQ(int(i/3), s.mol()[i]);
  mpart.push_back(21); mpart.push_back(22); mpart.push_back(23);

  s.delPart(mpart);
  mpart.clear();
  mpart = s.lastMolIDVec();
  EXPECT_EQ(294, mpart[0]);
  EXPECT_EQ(s.nMol(), int(s.mol2part().size())-1);
  EXPECT_EQ(s.natom(), s.mol2part().back());
  EXPECT_EQ(s.nMol() - 1, s.mol().back());
  EXPECT_EQ(s.nMol() - 1, s.mol().back());
  for (int i=0;i<s.natom();++i) EXPECT_EQ(int(i/3), s.mol()[i]);

  s.addMol("../forcefield/data.spce");
  mpart.clear();
  mpart = s.lastMolIDVec();
  EXPECT_EQ(297, mpart[0]);
  mpart.clear();
  mpart = s.lastMolIDVec();
  EXPECT_EQ(297, mpart[0]);
  EXPECT_EQ(s.nMol(), int(s.mol2part().size())-1);
  EXPECT_EQ(s.natom(), s.mol2part().back());
  EXPECT_EQ(s.nMol() - 1, s.mol().back());
  for (int i=0;i<s.natom();++i) EXPECT_EQ(int(i/3), s.mol()[i]);

  // delete all particles
  s.delPart(s.listAtoms());
  EXPECT_EQ(0, s.natom());
  EXPECT_EQ(0, s.nMol());
  EXPECT_EQ(s.nMol(), int(s.mol2part().size())-1);
  EXPECT_EQ(s.natom(), s.mol2part().back());

  s.addMolInit("../forcefield/data.atom");
  s.addMol("../forcefield/data.atom");
  EXPECT_EQ(1, s.natom());
  EXPECT_EQ(1, s.nMol());
  EXPECT_EQ(s.nMol(), int(s.mol2part().size())-1);
  EXPECT_EQ(s.natom(), s.mol2part().back());
  EXPECT_EQ(s.nMol() - 1, s.mol().back());
  for (int i=0;i<s.natom();++i) EXPECT_EQ(int(i/3), s.mol()[i]);
  EXPECT_EQ(1, s.checkSizes());
}

TEST(Space, xMolGen) {
  Space s(3);
  s.readXYZBulk(3, "water", "../unittest/spce/srsw/spce_sample_config_periodic1.xyz");
  s.xMolGen();
  EXPECT_EQ(100, s.nMol());
  for (int iMol = 0; iMol < s.nMol(); ++iMol) {
    EXPECT_EQ(3, int(s.xMol()[iMol].size()));
  }
  EXPECT_NEAR(-5.221309047080, s.xMol()[0][0][0], 1e-15);
  EXPECT_NEAR(-1.904156489730, s.xMol()[1][0][0], 1e-15);
  EXPECT_NEAR(-8.288503209960, s.xMol()[2][2][2], 1e-15);
  EXPECT_EQ(1, s.checkSizes());
}

TEST(Space, minl) {
  Space s(3);
  s.readXYZBulk(3, "water", "../unittest/spce/test52.xyz");
  EXPECT_NEAR(0, s.minl(), 1e-15);
  s.initBoxLength(5, 0);  s.initBoxLength(6, 1);  s.initBoxLength(7, 2);
  EXPECT_NEAR(5, s.minl(), 1e-15);
  EXPECT_EQ(1, s.checkSizes());
}

TEST(Space, nMol) {
  Space s(3);
  s.readXYZBulk(3, "water", "../unittest/spce/srsw/spce_sample_config_periodic1.xyz");
  s.addMolInit("../forcefield/data.spce");
  EXPECT_EQ(100, s.nMol());
  s.addMol("../forcefield/data.spce");
  EXPECT_EQ(101, s.nMol());
  vector<int> mpart(3);
  mpart[0] = 0; mpart[1] = 1; mpart[2] = 2;
  s.delPart(mpart);
  EXPECT_EQ(100, s.nMol());
  EXPECT_EQ(1, s.checkSizes());
}

TEST(Space, wrapMol) {
  Space s(3);
  const double boxl = 24.8586887;
  s.initBoxLength(boxl);
  s.readXYZBulk(3, "water", "../unittest/spce/test.xyz");
  ranInitByDate();
  vector<int> mpart = s.randMol();
  s.wrap(mpart);
  for (int dim = 0; dim < s.dimen(); ++dim) {
    for (int i = 0; i < int(mpart.size()); ++i) {
      s.xset(s.x(mpart[i],dim) + s.boxLength()[dim], mpart[i], dim);
    }
  }
  const double x1 = s.x(mpart[0],0);
  s.wrap(mpart);
  EXPECT_NEAR(boxl, fabs(x1 - s.x(mpart[0],0)), 1e-12);
  EXPECT_EQ(1, s.checkSizes());
}

//TEST(Space, avb) {
//  Space s(3);
//  const double boxl = 24.8586887;
//  for (int dim=0; dim < s.dimen(); ++dim) s.initBoxLength(boxl,dim);
//  s.readXYZBulk(3, "water", "../test/spce/test.xyz");
//  ranInitByDate();
//  const double rabove = 5, rbelow = 2.5;
//  for (int i = 0; i < 10; ++i) {
//    vector<int> mpart = s.randMol();
//    vector<int> mMol(1, s.mol()[mpart[0]]);
//    vector<int> jmpart = s.randMolDiff(mMol);
//    s.checkBond("spce", 1e-12);
//    int bonded = s.avb(mpart, jmpart, rabove, rbelow, "bonded");
//    s.checkBond("spce", 1e-12);
//    EXPECT_EQ(1, s.bonded(mpart, jmpart, rabove, rbelow));
//    bonded = s.avb(mpart, jmpart, rabove, rbelow, "nonbonded");
//    s.checkBond("spce", 1e-12);
//    EXPECT_EQ(1, bonded);
//    EXPECT_EQ(0, s.bonded(mpart, jmpart, rabove, rbelow));
//    const double xold = s.x(mpart[0],0);
//    s.xStore(mpart);
//    bonded = s.avb(mpart, jmpart, rabove, rbelow, "bonded");
//    s.checkBond("spce", 1e-12);
//    EXPECT_EQ(0, bonded);
//    s.restore(mpart);   // should restore to position before moving to bonded
//    EXPECT_NEAR(xold, s.x(mpart[0],0), 1e-12);
//  }
//  EXPECT_EQ(1, s.checkSizes());
//}

TEST(Space, randMolDiffANDSubset) {
  Space s(3);
  s.readXYZBulk(3, "water", "../unittest/spce/test52.xyz");
  ranInitByDate();
  for (int i = 0; i < 3*s.nMol(); ++i) {
    vector<int> mpart = s.randMol();
    vector<int> mMol(1, s.mol()[mpart[0]]);
    vector<int> jmpart = s.randMolDiff(mMol);
    EXPECT_NE(s.mol()[mpart[0]], s.mol()[jmpart[0]]);
    jmpart.clear();
    vector<int> jmMol(1, s.mol()[mpart[0]]);
    jmpart = s.randMolSubset(jmMol);
    EXPECT_EQ(s.mol()[mpart[0]], s.mol()[jmpart[0]]);
  }
  // build mpart as every molecule except one, then see if you select the only one
  vector<int> mMol;
  for (int i = 0; i < s.nMol(); ++i) {
    if (i != 4) mMol.push_back(i);
  }
  for (int i = 0; i < s.nMol(); ++i) {
    vector<int> jmpart = s.randMolDiff(mMol);
    EXPECT_EQ(4, s.mol()[jmpart[0]]);
  }
  EXPECT_EQ(1, s.checkSizes());
}

TEST(Space, listAtomsandMols) {
  Space s(3);
  s.readXYZBulk(3, "water", "../unittest/spce/test52.xyz");
  s.addMolInit("../forcefield/data.spce");
  EXPECT_EQ(1, s.checkSizes());
  s.addMol("../forcefield/data.spce");
  EXPECT_EQ(1, s.checkSizes());
  vector<int> mpart;
  mpart.push_back(0); mpart.push_back(1); mpart.push_back(2);
  s.delPart(mpart);
  EXPECT_EQ(1, s.checkSizes());
}

TEST(Space, cellListWater) {
  Space s(3);
  s.readXYZBulk(3, "water", "../unittest/spce/test52.xyz");
  EXPECT_EQ(0, s.cellType());
  s.updateCells(2,14);
  EXPECT_EQ(0, s.cellType());
  const double boxl = 24.8586887;
  s.initBoxLength(boxl);
  s.updateCells(2,boxl/2);
  EXPECT_EQ(0, s.cellType());
  const int ncell = 4;
  s.updateCells(boxl/double(ncell), 0.001);
  EXPECT_EQ(1, s.cellType());
  for (int dim=0; dim < s.dimen(); ++dim) EXPECT_EQ(ncell, s.nCellVec()[dim]);
  EXPECT_EQ(1, s.checkSizes());
}

// randomly update position of molecules with quaternions, and check that they did not change
TEST(Space, quat2posANDqMolInit) {
  Space s(3);
  s.readXYZBulk(3, "water", "../unittest/spce/test52.xyz");
  ranInitByDate();
  for (int i = 0; i < 10; ++i) {
    vector<int> mpart = s.randMol();
    const int iMol = s.mol()[mpart[0]];
    vector<vector<double> > xmol = s.xMol()[iMol];
    s.quat2pos(iMol);
    for (int i = 0; i < int(xmol.size()); ++i) {
      for (int dim = 0; dim < s.dimen(); ++dim) {
        EXPECT_NEAR(s.x(s.mol2part()[iMol]+i, dim), xmol[i][dim], 1e-16);
      }
    }
  }
  EXPECT_EQ(1, s.checkSizes());
}

TEST(Space, patch6) {
  Space s(3);
  s.initBoxLength(10);
  s.readXYZBulk(2, "onePatch", "../unittest/patch/onePatch6.xyz");
  EXPECT_EQ(s.nMol(), 6);
  EXPECT_NEAR(s.xMolRef()[0][0][0], 0, 1e-15);
  EXPECT_NEAR(s.xMolRef()[0][1][0], 1, 1e-15);
  EXPECT_NEAR(s.xMolRef()[1][0][0], 0, 1e-15);
  EXPECT_NEAR(s.xMolRef()[1][1][0], 1, 1e-15);
  s.settype(8, 2);
  s.settype(10, 2);
  EXPECT_EQ(0, s.type()[0]);
  EXPECT_EQ(1, s.type()[1]);
  EXPECT_EQ(2, s.type()[8]);
  EXPECT_EQ(3, int(s.nType().size()));
  EXPECT_EQ(4, s.nType()[0]);
  EXPECT_EQ(6, s.nType()[1]);
  EXPECT_EQ(2, s.nType()[2]);
  vector<int> mpart; mpart.push_back(0); mpart.push_back(1);
  s.delPart(mpart);
  EXPECT_EQ(3, int(s.nType().size()));
  EXPECT_EQ(3, s.nType()[0]);
  EXPECT_EQ(5, s.nType()[1]);
  EXPECT_EQ(2, s.nType()[2]);
  EXPECT_EQ(1, s.checkSizes());
}

TEST(Space, cellListOnePatch) {
  for (int cell = 0; cell < 2; ++cell) {
    ranInitByDate();
    Space s(3);
    s.readXYZBulk(2, "onePatch", "../unittest/patch/onePatch50.xyz");
    s.addMolInit("../forcefield/data.onePatch");
    EXPECT_EQ(0, s.cellType());
    const double rCut = 1.5;
    s.initBoxLength(10);
    s.initAtomCut(cell);
    s.updateCells(rCut, rCut);
    EXPECT_EQ(1, s.cellType());
    EXPECT_EQ(s.nCell(), int(s.neighCell().size()));
    for (int mi = 0; mi < s.nCell(); ++mi) {
      EXPECT_EQ(27, int(s.neighCell()[mi].size()));
    }
    EXPECT_EQ(6, s.nCellVec()[0]);
    EXPECT_EQ(6*6*6, s.nCell());
    if (cell == 0) {
      s.buildNeighListCell(0);
    } else {
      s.buildNeighListCellAtomCut(0);
    }
    //EXPECT_EQ(1, s.checkCellList());

    // add and remove particles
    vector<int> mpart = s.randMol();
    s.delPart(mpart);
  //  cout << "deling mol " << s.mol()[mpart.front()] << " mp " << mpart.front() << endl;
    s.addMol("../forcefield/data.onePatch");

    EXPECT_EQ(1, s.checkCellList());

    s.cellOff();
    EXPECT_EQ(0, s.cellType());
    //EXPECT_EQ(1, s.checkSizes());
  }
}

TEST(Space, readDataSPCE) {
  Space s(3);
  s.addMolInit("../forcefield/data.spce");
  vector<double> xAdd(s.dimen());
  s.xAdd = xAdd;
  s.addMol("../forcefield/data.spce");
  //s.initData("../forcefield/data.spce");
  EXPECT_EQ(3, s.natom());
  //EXPECT_NEAR(0.816496580927726000, s.x(1,0), 1e-15);
  EXPECT_NEAR(1, s.x(1,0), 1e-15);
  EXPECT_EQ(0, s.mol().at(2));
  EXPECT_EQ(1, s.type(2));
  EXPECT_EQ(s.natom()/3, s.nType()[0]);
  EXPECT_EQ(s.natom()*2/3, s.nType()[1]);
  EXPECT_EQ(2, s.nParticleTypes());
  EXPECT_EQ(0, s.moltype().front().compare("../forcefield/data.spce"));
  EXPECT_EQ(0, s.moltype().back().compare("../forcefield/data.spce"));
  EXPECT_EQ(1, s.checkSizes());
  shared_ptr<Space> stmp = s.findAddMolInList("../forcefield/data.spce");
  EXPECT_EQ(2, int(stmp->bondList().size()));
  EXPECT_EQ(0, stmp->bondList()[0][0]);
  EXPECT_EQ(0, stmp->bondList()[0][1]);
  EXPECT_EQ(1, stmp->bondList()[0][2]);
  EXPECT_EQ(0, stmp->bondList()[1][0]);
  EXPECT_EQ(0, stmp->bondList()[1][1]);
  EXPECT_EQ(2, stmp->bondList()[1][2]);
  EXPECT_EQ(1, int(stmp->bondParam().size()));
  EXPECT_NEAR(450, stmp->bondParam()[0][0], DTOL);
  EXPECT_NEAR(1, stmp->bondParam()[0][1], DTOL);
  vector<double> params = stmp->bondParams(0, 1);
  EXPECT_NEAR(450, params[0], DTOL);
  EXPECT_NEAR(1, params[1], DTOL);
  params = stmp->bondParams(1, 0);
  EXPECT_NEAR(450, params[0], DTOL);
  EXPECT_NEAR(1, params[1], DTOL);
  vector<vector<int> > bl = s.listBonds(0);
  EXPECT_EQ(2, int(bl.size()));
  EXPECT_EQ(0, bl[0][0]);
  EXPECT_EQ(0, bl[0][1]);
  EXPECT_EQ(1, bl[0][2]);
  EXPECT_EQ(0, bl[1][1]);
  EXPECT_EQ(2, bl[1][2]);
  bl = s.listBonds(1);
  EXPECT_EQ(1, int(bl.size()));
  bl = s.listBonds(2);
  EXPECT_EQ(1, int(bl.size()));

  EXPECT_EQ(1, int(stmp->angleList().size()));
  EXPECT_EQ(0, stmp->angleList()[0][0]);
  EXPECT_EQ(1, stmp->angleList()[0][1]);
  EXPECT_EQ(0, stmp->angleList()[0][2]);
  EXPECT_EQ(2, stmp->angleList()[0][3]);
  EXPECT_EQ(1, int(stmp->angleParam().size()));
  EXPECT_NEAR(55, stmp->angleParam()[0][0], DTOL);
  EXPECT_NEAR(109.47/180*PI, stmp->angleParam()[0][1], DTOL);
  params = stmp->angleParams(0, 1, 2);
  EXPECT_NEAR(55, params[0], DTOL);
  EXPECT_NEAR(109.47/180*PI, params[1], DTOL);
  params = stmp->angleParams(1, 0, 2);
  EXPECT_NEAR(55, params[0], DTOL);
  EXPECT_NEAR(109.47/180*PI, params[1], DTOL);
  params = stmp->angleParams(2, 0, 1);
  EXPECT_NEAR(55, params[0], DTOL);
  EXPECT_NEAR(109.47/180*PI, params[1], DTOL);
  params = stmp->angleParams(2, 1, 0);
  EXPECT_NEAR(55, params[0], DTOL);
  EXPECT_NEAR(109.47/180*PI, params[1], DTOL);
  vector<vector<int> > al = s.listAngles(0,1);
  EXPECT_EQ(1, int(al.size()));
  EXPECT_EQ(0, al[0][0]);
  EXPECT_EQ(1, al[0][1]);
  EXPECT_EQ(0, al[0][2]);
  EXPECT_EQ(2, al[0][3]);
  al = s.listAngles(0,2); EXPECT_EQ(1, int(al.size()));
  al = s.listAngles(1,2); EXPECT_EQ(1, int(al.size()));
  al = s.listAngles(2,1); EXPECT_EQ(1, int(al.size()));
  al = s.listAngles(1,0); EXPECT_EQ(1, int(al.size()));
}

TEST(Space, readDataEqultl43) {
int ntest = 1;
#ifdef JSON_
  ntest = 2;
#endif  // JSON_
  for (int itest = 0; itest < ntest; ++itest) {
    Space s(3);
    if (itest == 0) {
      s.initData("../forcefield/data.equltl43");
    } else {
      s.initData("../forcefield/equltl43.json");
    }
    EXPECT_EQ(4, s.natom());
    EXPECT_NEAR(4./3./sqrt(3), s.x(1,1), 1e-15);
    EXPECT_NEAR(4./3./2., s.x(2,0), 1e-15);
    EXPECT_NEAR(-4./3./2./sqrt(3), s.x(2,1), 1e-15);
    EXPECT_NEAR(-s.x(2,0), s.x(3,0), 1e-15);
    EXPECT_NEAR(s.x(2,1), s.x(3,1), 1e-15);
    EXPECT_NEAR(0, s.x(1,0), 1e-15);
    EXPECT_NEAR(0, s.x(1,2), 1e-15);
    EXPECT_NEAR(0, s.x(2,2), 1e-15);
    EXPECT_NEAR(0, s.x(3,2), 1e-15);
    EXPECT_EQ(0, s.mol().at(2));
    EXPECT_EQ(0, s.type().at(0));
    EXPECT_EQ(1, s.type().at(1));
    EXPECT_EQ(1, s.nType()[0]);
    EXPECT_EQ(3, s.nType()[1]);
    EXPECT_EQ(2, s.nParticleTypes());
    if (itest == 0) {
      EXPECT_EQ(0, s.moltype().front().compare("../forcefield/data.equltl43"));
      EXPECT_EQ(0, s.moltype().back().compare("../forcefield/data.equltl43"));
    } else {
      EXPECT_EQ(0, s.moltype().front().compare("../forcefield/equltl43.json"));
      EXPECT_EQ(0, s.moltype().back().compare("../forcefield/equltl43.json"));
    }
    EXPECT_EQ(1, s.checkSizes());
  }
}

TEST(Space, twoequltl43) {
  Space s(3);
  s.initBoxLength(9);
  s.readXYZBulk(4, "equltl43", "../unittest/equltl43/two.xyz");
  EXPECT_EQ(8, s.natom());
  EXPECT_NEAR(4./3./sqrt(3), s.x(1,1), 1e-15);
  EXPECT_NEAR(4./3./2., s.x(2,0), 1e-15);
  EXPECT_NEAR(-4./3./2./sqrt(3), s.x(2,1), 1e-15);
  EXPECT_NEAR(-s.x(2,0), s.x(3,0), 1e-15);
  EXPECT_NEAR(s.x(2,1), s.x(3,1), 1e-15);
  EXPECT_NEAR(0, s.x(1,0), 1e-15);
  EXPECT_NEAR(0, s.x(1,2), 1e-15);
  EXPECT_NEAR(0, s.x(2,2), 1e-15);
  EXPECT_NEAR(0, s.x(3,2), 1e-15);
  EXPECT_EQ(0, s.mol().at(2));
  EXPECT_EQ(0, s.type().at(0));
  EXPECT_EQ(1, s.type().at(1));
  EXPECT_EQ(2, s.nType()[0]);
  EXPECT_EQ(6, s.nType()[1]);
  EXPECT_EQ(2, s.nParticleTypes());
  EXPECT_EQ(0, s.moltype().front().compare("equltl43"));
  EXPECT_EQ(0, s.moltype().back().compare("equltl43"));
  EXPECT_EQ(1, s.checkSizes());
}

TEST(Space, clone) {
  Space s(3);
  s.initBoxLength(9);
  s.readXYZBulk(4, "equltl43", "../unittest/equltl43/two.xyz");
  shared_ptr<Space> s2 = s.cloneShrPtr();
}

TEST(Space, scaleMol) {
  Space s(3);
  s.readXYZBulk(4, "equltl43", "../unittest/equltl43/two.xyz");
  EXPECT_NEAR(s.x(1, 1), 0.769800358919501000, 1e-15);
  EXPECT_NEAR(s.x(2,0), 0.666666666666666666, 1e-15);
  EXPECT_NEAR(s.x(2,1), -0.384900179459751000, 1e-15);
  vector<double> bondLength(4, 2./3./sqrt(3));
  s.scaleMol(0, bondLength);
  EXPECT_NEAR(s.x(1,1), 0.769800358919501000/2., 1e-15);
  EXPECT_NEAR(s.x(2,0), 0.666666666666666666/2., 1e-15);
  EXPECT_NEAR(s.x(2,1), -0.384900179459751000/2., 1e-15);
}

TEST(Space, writeRestart) {
  double xold, xold2;
  { Space s(3);
    s.initBoxLength(9);
    s.addMolInit("../forcefield/data.equltl43");
    s.addMol("../forcefield/data.equltl43");
    s.addMol("../forcefield/data.equltl43");
    s.tagAtom(2);
    s.initPBC(0, 0);
    s.writeRestart("tmp/rst");
    xold = s.x(0,0);
    xold2 = s.x(2,2);
  }

  Space s("tmp/rst");
  EXPECT_EQ(2, s.tag()[0]);
  EXPECT_EQ(8, s.natom());
  EXPECT_NEAR(xold, s.x(0,0), 1e-16);
  EXPECT_NEAR(xold2, s.x(2,2), 1e-16);
  EXPECT_EQ(s.periodic(0), 0);
}

TEST(Space, readxyzmulti) {
  Space s(3);
  int iConf = 0;
  s.readXYZ("../unittest/xyz1");
  while (!s.xyzFileEOF()) {
  //while (s.xyzFile()->eof() != std::istream::traits_type::eof()) {
    EXPECT_EQ(12+iConf, s.natom());
    EXPECT_NEAR(0.01*iConf, s.x(0,0), 1e-15);
    s.readXYZ("../unittest/xyz1");
    ++iConf;
  }
  EXPECT_EQ(2, iConf);
}

TEST(Space, delTypePart) {
  Space s(3);
  s.readXYZ("../unittest/cg3_60_43_1snap.xyz");
  s.initBoxLength(12);
  EXPECT_EQ(s.natom(), 900);
  EXPECT_EQ(s.nParticleTypes(), 3);
  s.delTypePart(2);
  EXPECT_EQ(s.natom(), 450);
  EXPECT_EQ(s.nParticleTypes(), 3);
  s.delTypePart(0);
  EXPECT_EQ(s.natom(), 225);
  EXPECT_EQ(s.nParticleTypes(), 3);
}

TEST(Space, floodFill) {
  Space s(3);
  s.addMolInit("../forcefield/data.cg3_60_43_1");
  std::ifstream file("../unittest/cg3_60_43_1snap.xyz");
  s.readXYZ(file);
  s.addTypeForCluster(1);
  s.updateClusters(4./3.);
  EXPECT_EQ(12, s.nClusters());
}

TEST(Space, swapPositionsANDstoreAll) {
  Space s(3);
  s.addMolInit("../forcefield/data.cg3_60_43_1");
  for (int i = 0; i < 10; ++i) s.addMol("../forcefield/data.cg3_60_43_1");
  s.initBoxLength(12);
  shared_ptr<Space> s2 = s.cloneShrPtr();
  s2->xStoreAll();
  const double posold = s.x(2, 2);
  vector<int> mpart;
  for (int i = 0; i < 4; ++i) mpart.push_back(i);
  s.randDisp(mpart, 2);
  s.randRotate(mpart, 2);
  const double posnew = s.x(2, 2);
  EXPECT_NE(posnew, posold);
  s.swapPositions(s2.get());
  EXPECT_NEAR(posnew, s2->x(2, 2), 1e-14);
  s2->restoreAll();
  EXPECT_NEAR(posold, s2->x(2, 2), 1e-14);
}

TEST(Space, maxMolDist) {
  Space s(3);
  s.initBoxLength(12);
  s.addMolInit("../forcefield/data.cg3_60_43_1");
  EXPECT_NEAR(0.7698, s.maxMolDist(), 1e-4);
  s.addMol("../forcefield/data.cg3_60_43_1");
  s.addMol("../forcefield/data.cg3_60_43_1");
  s.addMol("../forcefield/data.cg3_60_43_1");
  s.addMol("../forcefield/data.cg3_60_43_1");
  s.addMol("../forcefield/data.cg3_60_43_1");
  EXPECT_NEAR(0.7698, s.maxMolDist(), 1e-4);
}

TEST(Space, readXYZ) {
  Space s(3);
  s.addMolInit("../forcefield/data.cg3_60_43_1");
  std::ifstream file("../unittest/cg3_60_43_1snap.xyz");
  s.readXYZ(file);
  EXPECT_EQ(900, s.natom());
  s.readXYZ(file);
  EXPECT_EQ(900, s.natom());
  //EXPECT_EQ(900+896, s.natom());
  s.readXYZ(file);
  EXPECT_TRUE(file.eof());
}

TEST(Space, imol2mpart) {
  Space s(3);
  s.addMolInit("../forcefield/data.equltl43");
  s.addMol("../forcefield/data.equltl43");
  s.addMol("../forcefield/data.equltl43");
  const vector<int> mpart = s.imol2mpart(0);
  EXPECT_EQ(4, int(mpart.size()));
  EXPECT_EQ(0, mpart[0]);
  EXPECT_EQ(3, mpart[3]);
  const double xold = s.x(0,0);
  const vector<double> r(3, 1.);
  s.transMol(0, r);
  EXPECT_NEAR(xold+1., s.x(0,0), 1e-14);
}

TEST(Space, randRotateMulti) {
  Space s(3);
  s.initBoxLength(12);
  s.addMolInit("../forcefield/data.cg3_60_43_1");
  std::ifstream file("../unittest/cg3_60_43_1snap.xyz");
  s.readXYZ(file);

  s.addTypeForCluster(1);
  s.updateClusters(4./3.);
  EXPECT_EQ(12, s.nClusters());

  // rigidly translate clusters
  //s.printXYZ("tm1234.xyz", 1);
  s.xClusterGen();
  for (int i = 2; i < 3; ++i) {
    s.randRotateMulti(s.clusterList()[i], 10);
  }
  //s.printXYZ("tm1234.xyz", 0);

// HWH eigenvalue  // compute shape of clusters
//  s.xClusterShape();
//  EXPECT_NEAR(s.clusterAsphericityAv(), 0.6204511597901633, 1e-5);
//  EXPECT_NEAR(s.clusterAcylindricityAv(), 0.37957508233844489, 1e-5);
//  EXPECT_NEAR(s.clusterRelShapeAnisoAv(), 0.026661790962870163, 1e-5);
//  EXPECT_NEAR(s.clusterRgAv(), 2.1327427724005301, 1e-5);
}

TEST(Space, inertialTensor) {
  Space s(3);
  s.initBoxLength(12);
  s.addMolInit("../forcefield/data.cg3_60_43_1");
  s.addMol("../forcefield/data.cg3_60_43_1");
  vector<int> mpart;
  vector<vector<double> > x(s.natom(), vector<double>(s.dimen()));
  for (int i = 0; i < s.natom(); ++i) {
    mpart.push_back(i);
  }
  vector<double> rcm = s.rcom(mpart);
  for (int i = 0; i < s.natom(); ++i) {
    for (int dim = 0; dim < s.dimen(); ++dim) {
      x[i][dim] = s.x(i, dim) - rcm[dim];
      //cout << " x "  << x[i][dim] << endl;
    }
  }
  vector<vector<double> > tensor = s.inertialTensor(mpart);
  vector<vector<double> > evectors;
  vector<double> evalues;
//HWH eignenvalue  jacobi(tensor, evalues, evectors);
//  EXPECT_NEAR(2*0.88888888888+1.777777777, evalues[0]+evalues[1]+evalues[2], 1e-9);
//  vector<vector<double> > xref = matMul(x, evectors);
}

#ifdef XDRFILE_H_
TEST(Space, readXTC) {
  Space s(3);
  s.addMolInit("../forcefield/data.lj");
  for (int i =0; i < 245; ++i) {
    s.addMol("../forcefield/data.lj");
  }
  XDRFILE* trjFileXDR;
  XDRFILE* trjFileXDR2;
  trjFileXDR = xdrfile_open("../unittest/1L2Y.xtc", "r");
  trjFileXDR2 = xdrfile_open("tmp/1L2Ycpy.xtc", "w");
  EXPECT_EQ(0, s.readXTC("../unittest/1L2Y.xtc", trjFileXDR));
  s.writeXTC(trjFileXDR2);
  EXPECT_NEAR(0.8275, s.x(0,0), 1e-5);
  EXPECT_EQ(1, s.readXTC("../unittest/1L2Y.xtc", trjFileXDR));
  s.writeXTC(trjFileXDR2);
  xdrfile_close(trjFileXDR);
  xdrfile_close(trjFileXDR2);
  XDRFILE* trjFileXDR3;
  trjFileXDR3 = xdrfile_open("tmp/1L2Ycpy.xtc", "r");
  EXPECT_EQ(0, s.readXTC("tmp/1L2Ycpy.xtc", trjFileXDR3));
  EXPECT_NEAR(0.8275, s.x(0,0), 1e-3);
  EXPECT_EQ(0, s.readXTC("tmp/1L2Ycpy.xtc", trjFileXDR3));
  EXPECT_NEAR(0., s.x(0,0), 1e-5);
  EXPECT_EQ(1, s.readXTC("tmp/1L2Ycpy.xtc", trjFileXDR3));
  xdrfile_close(trjFileXDR3);
}
#endif

TEST(Space, xAddANDqMolAlt) {
  Space s(3);
  s.initBoxLength(1e3);
  s.addMolInit("../forcefield/data.cg3_60_1_1");
  vector<double> xAdd(s.dimen());
  s.xAdd = xAdd;
  s.addMol("../forcefield/data.cg3_60_1_1");
  EXPECT_EQ(s.x(0,0), 0);
  EXPECT_EQ(s.x(0,1), 0);
  EXPECT_EQ(s.x(0,2), 0);
  EXPECT_EQ(s.x(1,0), 0);
  EXPECT_EQ(s.x(1,1), 0.577350269189626000);
  EXPECT_EQ(s.x(1,2), 0);
  EXPECT_EQ(s.x(2,0), 0.5);
  EXPECT_EQ(s.x(2,1), -0.288675134594813000);
  EXPECT_EQ(s.x(2,2), 0);
  EXPECT_EQ(s.x(3,0), -0.5);
  EXPECT_EQ(s.x(3,1), -0.288675134594813000);
  EXPECT_EQ(s.x(3,2), 0);
  xAdd[s.dimen()-1] = 1;
  s.xAdd = xAdd;
  s.addMol("../forcefield/data.cg3_60_1_1");
  EXPECT_EQ(s.x(4,0), 0);
  EXPECT_EQ(s.x(4,1), 0);
  EXPECT_EQ(s.x(4,2), 1);
  EXPECT_EQ(s.x(5,0), 0);
  EXPECT_EQ(s.x(5,1), 0.577350269189626000);
  EXPECT_EQ(s.x(5,2), 1);
  EXPECT_EQ(s.x(6,0), 0.5);
  EXPECT_EQ(s.x(6,1), -0.288675134594813000);
  EXPECT_EQ(s.x(6,2), 1);
  EXPECT_EQ(s.x(7,0), -0.5);
  EXPECT_EQ(s.x(7,1), -0.288675134594813000);
  EXPECT_EQ(s.x(7,2), 1);
  xAdd[0] = 0;
  xAdd[1] = 1;
  xAdd[2] = 0;
  s.xAdd = xAdd;
  s.addMol("../forcefield/data.cg3_60_1_1");
  s.qMolAlt(2, 0, 1);
  s.qMolAlt(2, 3, 0);
  s.quat2pos(2);
  EXPECT_EQ(s.x(8,0), 0);
  EXPECT_EQ(s.x(8,1), 1);
  EXPECT_EQ(s.x(8,2), 0);
  EXPECT_EQ(s.x(9,0), 0);
  EXPECT_EQ(s.x(9,1), 1-0.577350269189626000);
  EXPECT_EQ(s.x(9,2), 0);
  EXPECT_EQ(s.x(10,0), 0.5);
  EXPECT_EQ(s.x(10,1), 1+0.288675134594813000);
  EXPECT_EQ(s.x(10,2), 0);
  EXPECT_EQ(s.x(11,0), -0.5);
  EXPECT_EQ(s.x(11,1), 1+0.288675134594813000);
  EXPECT_EQ(s.x(11,2), 0);
}

TEST(Space, roundSquares) {
  Space s(2);
  s.initBoxLength(10);
  s.addMolInit("../forcefield/data.onePatch");
//  vector<double> xAdd(s.dimen());
//  s.xAdd = xAdd;
  s.addMol("../forcefield/data.onePatch");
  vector<int> mpart(1);
  s.randRotate(mpart, 100);

}

TEST(Space, setAtomInBranch) {
  Space s(3);
  s.addMolInit("../forcefield/data.lj");
  vector<double> xAdd(s.dimen()), r1(s.dimen()), r2(s.dimen()), r3(s.dimen());
  xAdd[0] = 0.5; xAdd[1] = 0.4; xAdd[2] = 0.768114575;
  r1 = xAdd;
  s.xAdd = xAdd;
  s.addMol("../forcefield/data.lj");
  xAdd[0] = -0.8; xAdd[1] = -0.2; xAdd[2] = 0.565685425;
  r2 = xAdd;
  s.xAdd = xAdd;
  s.addMol("../forcefield/data.lj");
  xAdd[0] = 0; xAdd[1] = 0; xAdd[2] = 0;
  s.xAdd = xAdd;
  s.addMol("../forcefield/data.lj");
  s.xAdd = xAdd;
  s.addMol("../forcefield/data.lj");
  s.setAtomInBranch(0, 1, 2, 3, 1.61630081, 1.61630081, 1);
  r3[0] = s.x(2, 0); r3[1] = s.x(2, 1); r3[2] = s.x(2,2);
  EXPECT_NEAR(cos(1.61630081), vecDotProd(r1, r3), 1e-8);
  EXPECT_NEAR(cos(1.61630081), vecDotProd(r1, r2), 1e-8);
  EXPECT_NEAR(cos(1.61630081), vecDotProd(r2, r3), 1e-8);
}

TEST(Space, randPosition) {
  Space s(3);
  s.initBoxLength(10);
  s.addMolInit("../forcefield/data.cg3_60_1_1");
  vector<double> x(s.dimen(), 0.);
  s.xAdd = x;
  s.addMol("../forcefield/data.cg3_60_1_1");
  for (int i = 0; i < 10; ++i) {
    const double maxMove = 0.1;
    const vector<double> r = s.randPosition(0, maxMove);
//    cout << "r " << r[0] << " " << r[1] << " " << r[2] << endl;
    const double max = *std::max_element(r.begin(), r.begin()+int(r.size()));
    const double min = *std::min_element(r.begin(), r.begin()+int(r.size()));
    EXPECT_LE(-maxMove, min);
    EXPECT_GE(maxMove, max);
  }
}

TEST(Space, scaleDomain) {
  Space s(3);
  s.initBoxLength(10);
  s.addMolInit("../forcefield/data.cg3_60_1_1");
  vector<double> x(s.dimen(), 0.);
  s.xAdd = x;
  s.addMol("../forcefield/data.cg3_60_1_1");
  s.addMol("../forcefield/data.cg3_60_1_1");
  s.scaleDomain(1.1);
  EXPECT_EQ(1, s.checkBond(10*DTOL));
}

//TEST(Space, modBondAngle) {
//  Space s2(3);
//  for (int i = 0; i < 100; ++i) {
//    Space s(3);
//    s.addMolInit("../forcefield/data.cg3_60_1_1");
//    vector<double> xAdd(s.dimen());
//    s.xAdd = xAdd;
//    s.addMol("../forcefield/data.cg3_60_1_1");
//    const vector<int> mpart = s.randMol();
//    s.randRotate(mpart, -1);
//    s.modBondAngle(0, PI/2, "../forcefield/data.cg3_60_1_1");
//  }
//}

TEST(Space, Q6) {

  // imperfect square lattice
  { Space sSqIm(2);
    sSqIm.addMolInit("../forcefield/data.onePatch1000");
    std::ifstream file("../unittest/roundSquare/square.xyz");
    sSqIm.readXYZ(file);
    EXPECT_NEAR(0.58613463576249036, sSqIm.Q6(1.), 1e-14);
  }

  // imperfect hex lattice
  { Space sHexIm(2);
    sHexIm.addMolInit("../forcefield/data.onePatch1000");
    std::ifstream file("../unittest/roundSquare/hex.xyz");
    sHexIm.readXYZ(file);
    EXPECT_NEAR(0.69295149889390972, sHexIm.Q6(1.), 1e-14);
  }

  // imperfect rhombic lattice
  { Space sRhmIm(2);
    sRhmIm.addMolInit("../forcefield/data.onePatch1000");
    std::ifstream file("../unittest/roundSquare/lambda1.xyz");
    sRhmIm.readXYZ(file);
    EXPECT_NEAR(0.71040413847700912, sRhmIm.Q6(1.), 1e-14);
  }

  // perfect cubic lattice
  { Space sCub(3);
    sCub.addMolInit("../forcefield/data.atom");
    std::ifstream file("../unittest/lattice/cubic.xyz");
    sCub.readXYZ(file);
    EXPECT_NEAR(0.3535533905932739, sCub.Q6(1.1), 1e-14);
  }

  // perfect square lattice
  { Space sSq(2);
    sSq.addMolInit("../forcefield/data.atom");
    std::ifstream file("../unittest/lattice/square.xyz");
    sSq.readXYZ(file);
    EXPECT_NEAR(0.58630196997792861, sSq.Q6(1.1), 1e-14);
  }

  // perfect hex lattice
  { Space sHex(2);
    sHex.addMolInit("../forcefield/data.atom");
    std::ifstream file("../unittest/lattice/hex.xyz");
    sHex.readXYZ(file);
    EXPECT_NEAR(0.74082934944560608, sHex.Q6(1.1), 1e-14);
  }

  // perfect fcc lattice
  { Space sFCC(3);
    sFCC.addMolInit("../forcefield/data.atom");
    std::ifstream file("../unittest/lattice/fcc.xyz");
    sFCC.readXYZ(file);
    EXPECT_NEAR(0.5745242597140704, sFCC.Q6(sqrt(2)/2+0.001), 1e-14);
  }
}

TEST(Space, xyTilt) {

  Space s(2);
  s.initBoxLength(12);
  s.setXYTilt(12);
  s.addMolInit("../forcefield/data.onePatch1000");
  for (int i = 0; i < 1000; ++i) s.addMol("../forcefield/data.onePatch1000");
  //s.printXYZ("tilt",1);
  EXPECT_NEAR(s.minBondLength(), 1, 10*DTOL);
}

TEST(Space, SQ) {
  Space s(3);
  s.initBoxLength(12);
  s.addMolInit("../forcefield/data.cg7mab1");
  EXPECT_NEAR(s.minBondLength(), 1.8122340298918349, 10*DTOL);
}

//TEST(Space, pos2euler) {
//  Space s(3);
//  for (int dim=0; dim < s.dimen(); ++dim) s.initBoxLength(12,dim);
//  s.addMolInit("../forcefield/data.twoPatch");
//  s.addMol();
//  shared_ptr<Space> s2 = s.cloneShrPtr();
//  s2->pos2euler(0);
//  for (int iatom = 0; iatom < s.natom(); ++iatom) {
//    for (int dim = 0; dim < s.dimen(); ++dim) {
//      EXPECT_NEAR(s.x(iatom,dim), s2->x(iatom,dim), DTOL);
//    }
//  }
//}

TEST(Space, randMolofType) {
  Space s(3);
  s.addMolInit("../forcefield/data.lj");
  s.addMolInit("../forcefield/data.ljb");
  s.addMol("../forcefield/data.lj");
  s.addMol("../forcefield/data.ljb");
  s.addMol("../forcefield/data.lj");
  s.addMol("../forcefield/data.ljb");
  int n=200, count = 0;
  for (int i = 0; i < n; ++i) {
    int iMol = s.randMolofType(0);
    EXPECT_TRUE( (iMol == 0) || (iMol == 2) );
    if (iMol == 0) ++count;
    iMol = s.randMolofType(1);
    EXPECT_TRUE( (iMol == 1) || (iMol == 3) );
  }
  EXPECT_NEAR(double(count)/double(n), 0.5, 0.15);
}

TEST(Space, molid) {
  Space s(3);
  EXPECT_EQ(int(s.molid().size()), 0);
  EXPECT_EQ(int(s.nMolType().size()), 1);
  s.addMolInit("../forcefield/data.lj");
  s.addMolInit("../forcefield/data.ljb");

  // add first molecule
  s.addMol("../forcefield/data.lj");
  EXPECT_EQ(int(s.molid().size()), 1);
  EXPECT_EQ(int(s.nMolType().size()), 2);
  EXPECT_EQ(s.molid()[0], 0);
  EXPECT_EQ(s.nMolType()[0], 1);
  EXPECT_EQ(s.nMolType()[1], 0);

  // add second molecule
  s.addMol("../forcefield/data.ljb");
  EXPECT_EQ(int(s.molid().size()), 2);
  EXPECT_EQ(int(s.nMolType().size()), 2);
  EXPECT_EQ(s.molid()[0], 0);
  EXPECT_EQ(s.molid()[1], 1);
  EXPECT_EQ(s.nMolType()[0], 1);
  EXPECT_EQ(s.nMolType()[1], 1);

  s.addMol("../forcefield/data.lj");
  s.addMol("../forcefield/data.ljb");
}

TEST(Space, readDataCG7MabAniso) {
int ntest = 1;
#ifdef JSON_
  ntest = 2;
#endif  // JSON_
  for (int itest = 0; itest < ntest; ++itest) {
//    cout << "itest " << itest << endl;
    Space s(3);
    if (itest == 0) {
      s.initData("../forcefield/data.cg7mabaniso");
    } else {
      s.initData("../forcefield/cg7mabAniso.json");
    }
    EXPECT_EQ(7, s.natom());
    EXPECT_NEAR(s.x(0,0), 0, 1e-15);
    EXPECT_NEAR(s.x(1,1), 1.3, 1e-15);
    EXPECT_NEAR(s.x(2,1), 1.3, 1e-15);
    EXPECT_NEAR(s.x(3,0), -1.060660172, 1e-15);
    EXPECT_NEAR(s.x(4,1), -1.767766953, 1e-15);
    EXPECT_NEAR(s.x(5,2), 0, 1e-15);
    EXPECT_NEAR(s.x(6,0), -s.x(4,1), 1e-15);
    EXPECT_NEAR(0, s.x(1,0), 1e-15);
    EXPECT_NEAR(0, s.x(1,2), 1e-15);
    EXPECT_NEAR(1, s.x(2,2), 1e-15);
    EXPECT_NEAR(0, s.x(3,2), 1e-15);
    EXPECT_EQ(0, s.mol().at(2));
    EXPECT_EQ(0, s.type().at(0));
    EXPECT_EQ(1, s.type().at(1));
    EXPECT_EQ(4, s.nType()[0]);
    EXPECT_EQ(1, s.nType()[1]);
    EXPECT_EQ(2, s.nType()[2]);
    EXPECT_EQ(3, s.nParticleTypes());
    for (int ibond = 0; ibond < 4; ++ibond) {
      EXPECT_EQ(-1, s.bondParam()[ibond][0]);
    }
    EXPECT_EQ(1.3, s.bondParam()[0][1]);
    EXPECT_EQ(1.5, s.bondParam()[2][1]);
    EXPECT_EQ(1.0, s.bondParam()[3][1]);
    EXPECT_EQ(3, s.bondList()[1][2]);
    EXPECT_EQ(135./180.*PI, s.angleParam()[0][1]);
    EXPECT_EQ(-1, s.angleParam()[2][0]);
    EXPECT_EQ(90./180.*PI, s.angleParam()[3][1]);
    EXPECT_EQ(5, s.angleList()[1][3]);
    EXPECT_EQ(6, s.angleList()[5][3]);
    EXPECT_EQ(5, s.angleList()[5][0]);
    if (itest == 0) {
      EXPECT_EQ(0, s.moltype().front().compare("../forcefield/data.cg7mabaniso"));
      EXPECT_EQ(0, s.moltype().back().compare("../forcefield/data.cg7mabaniso"));
    } else {
      EXPECT_EQ(0, s.moltype().front().compare("../forcefield/cg7mabAniso.json"));
      EXPECT_EQ(0, s.moltype().back().compare("../forcefield/cg7mabAniso.json"));
    }
    EXPECT_EQ(1, s.checkSizes());
  }
}

//// visualize this test
//TEST(Space, randRotateMaxDisp) {
//  ranInitByDate();
//  Space space(3);
//  space.addMolInit("../forcefield/data.onePatch");
//  vector<double> xAdd(3, 0.);
////  space.xAdd = xAdd;
//  space.addMol("../forcefield/data.onePatch");
//  vector<int> imove(2);
//  imove[0] = 0; imove[1] = 1;
//  for (int trial = 0; trial < 100; ++trial) {
//  //for (int trial = 0; trial < 10000; ++trial) {
//    space.xStore(imove);
//    space.randRotate(imove, 0.9);
//    cout << space.x(1, 0) << " " << space.x(1, 1) << " " << space.x(1, 2) << endl;
//    space.restore(imove);
//  }
//}

TEST(Space, replicate) {
  Space space(3);
  space.initBoxLength(9);
  space.addMolInit("../forcefield/data.cg3_180_43");
  for (int i = 0; i < 20; ++i) {
    space.addMol("../forcefield/data.cg3_180_43");
  }

  const double rClusterCut = 1.5;
  //const double rClusterCut = 2.5;

  space.addTypeForCluster(0);
  space.updateClusters(rClusterCut);
//  space.printClusterStat("hi");

//  space.printXYZ("hi", 1);
  shared_ptr<Space> spaceBig = space.cloneShrPtr();
  spaceBig->replicate();
//  spaceBig->printXYZ("hibig", 1);
  EXPECT_EQ(pow(2, space.dimen())*space.natom(), spaceBig->natom());

  spaceBig->updateClusters(rClusterCut);
//  spaceBig->printClusterStat("hibig");

  // for small rClusterCut, check that the number of clusters increases
  // by a factor of 2^3 for replication, which means that the clusters
  // do not percolate
  EXPECT_EQ(pow(2, space.dimen())*space.nClusters(), spaceBig->nClusters());
}

TEST(Space, argparse) {
  {
    Space space(3, {{"id", "123"}});
    EXPECT_EQ(space.dimen(), 3);
    EXPECT_EQ(space.id(), 123);
  }

  {
    Space space(3, {{"boxLength", "8.5"}});
    EXPECT_EQ(space.dimen(), 3);
    EXPECT_NEAR(space.minl(), 8.5, feasst::DTOL);
  }

  try {
    feasst::Space space(3, {{"not/a/proper/arg", "error"}});
    CATCH_PHRASE("is not recognized");
  }

  {
    auto space = makeSpace(2, {{"boxLength", "0.5"}});
    EXPECT_EQ(2, space->dimen());
    EXPECT_EQ(0.5, space->minl());
  }

  {
    auto space = makeSpace({{"boxLength", "8.5"}});
    EXPECT_EQ(3, space->dimen());
    EXPECT_EQ(8.5, space->minl());
  }
}

TEST(Space, useMaxRotAsAngle) {
  feasst::ranInitByDate();
  auto space = feasst::makeSpace({{"dimen", "3"}});
  space->addMolInit("../forcefield/data.spce");
  vector<double> xAdd(space->dimen());
  space->xAdd = xAdd;
  space->addMol();
  EXPECT_NEAR(space->x(2, 0), -0.333313247568237000, feasst::DTOL);
  EXPECT_NEAR(space->x(2, 1), 0.942816142731718000, feasst::DTOL);
  EXPECT_NEAR(space->x(2, 2), 0., feasst::DTOL);
  vector<double> axis = {1, 0, 0};
  space->randRotateAboutAxis(space->imol2mpart(0), axis, feasst::PI/2., 1);
  EXPECT_NEAR(space->x(2, 0), -0.333313247568237000, feasst::DTOL);
  EXPECT_NEAR(space->x(2, 1), 0., feasst::DTOL);
  EXPECT_NEAR(space->x(2, 2), -0.942816142731718000, feasst::DTOL);
  space->randRotateAboutAxis(space->imol2mpart(0), axis, -feasst::PI/2.);
  EXPECT_NE(space->x(2, 1), 0.942816142731718000);
}

TEST(Space, pointReflect) {
  Space s(3);
  s.initBoxLength(10);
  s.addMolInit("../forcefield/data.cg3_60_1_1");
  s.addMol("../forcefield/data.cg3_60_1_1");
  vector<double> x(s.dimen(), 0.);
  s.addMolInit("../forcefield/data.atom");
  s.xAdd = x;
  s.addMol("../forcefield/data.atom");
  //const vector<double> x = s.randPosition();
  s.printXYZ("ttt",1);
  s.pointReflectMol(0, x);
  s.printXYZ("ttt",0);
}

TEST(Space, lineReflect) {
  Space s(3);
  s.initBoxLength(20);
  s.addMolInit("../forcefield/data.cg3_60_1_1");
  vector<double> x1(s.dimen(), 0.), x2;
  x1[0] = 5.;
  s.xAdd = x1;
  s.addMol("../forcefield/data.cg3_60_1_1");
  s.addMolInit("../forcefield/data.atom");
  x1[0] = 0;
  s.xAdd = x1;
  s.addMol("../forcefield/data.atom");
  x2 = x1;
  x2[1] = 1;
  s.xAdd = x2;
  s.addMol("../forcefield/data.atom");
  //const vector<double> x = s.randPosition();
  s.printXYZ("ttt",1);
  s.lineReflectMol(0, x1, x2);
  s.printXYZ("ttt",0);
  s.pointReflectMol(0, x1);
  s.printXYZ("ttt",0);
  s.randLineReflectMolByBond(0);
  s.printXYZ("ttt",0);
}
