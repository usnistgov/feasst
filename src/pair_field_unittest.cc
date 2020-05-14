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
#include "pair_field_slit_lj.h"
#include "pair_field_slit_sw.h"
#include "pair_field_slit_linear.h"
#include "pair_field_sphere_sw.h"
#include "pair_field_cylinder_sw.h"
#include "pair_hybrid.h"
#include "mc.h"

void translateAndExpect(const vector<double> xPosition, const double peExpected,
  std::shared_ptr<feasst::Pair> pair) {
  std::vector<double> xAdd(pair->space()->dimen());
  for (int dim = 0; dim < pair->space()->dimen(); ++dim) {
    xAdd[dim] = xPosition[dim] - pair->space()->x(0, dim);
  }
  pair->space()->transMol(0, xAdd);
  pair->initEnergy();
  EXPECT_NEAR(peExpected, pair->peTot(), feasst::DTOL);
}

void checkRestart(std::shared_ptr<feasst::Pair> pair, const std::string des) {
  const double pe = pair->peTot();
  std::stringstream ss;
  ss << "tmp/fieldrst" << des;
  pair->writeRestart(ss.str().c_str());
  auto pair3 = makePair(pair->space(), ss.str().c_str());
  pair3->initEnergy();
  cout << "checkig rest" << endl;
  EXPECT_NEAR(pe, pair3->peTot(), feasst::DTOL);
  cout << "checkig rest" << endl;
}


void addMols(std::shared_ptr<feasst::Pair> pair, const int nMol) {
  auto criteria = feasst::makeCriteriaMetropolis({{"beta", "0.1"}});
  feasst::MC mc(pair, criteria);
  mc.nMolSeek(300);
}

TEST(PairFieldSlitLJ, lj) {
  auto space = feasst::makeSpace({{"dimen", "2"}, {"boxLength", "14"}});
  auto pair1 = feasst::makePairFieldSlitLJ(space);
  pair1->initData("../forcefield/data.lj");
  pair1->initSlit(0, 5, -5);
  auto pair2 = pair1->clone(space.get());
  pair1->initAlpha(9);
  pair1->initEps(0, 2./15.);
  pair1->initSig(0, 1.25);
  pair1->initDelta(0.25);
  pair2->initAlpha(3);
  pair2->initEps(0, -1.);
  pair2->initSig(0, 1.2);
  pair2->initDelta(0.15);
  auto pair = feasst::makePairHybrid(space.get());
  pair->addPair(pair1.get());
  pair->addPair(pair2);

  std::vector<double> xAdd(space->dimen(), 0.);
  pair->addMol(xAdd);
//  const double peExpect = 2*(2./15*pow(1.25/(5), 9) - pow(1.2/(5), 3));
  const double peExpect = 2*(2./15*pow(1.25/(5 + 0.25), 9) - pow(1.2/(5 + 0.15), 3));

  xAdd[0] = 7;
  translateAndExpect(xAdd, NUM_INF*2, pair);
  xAdd[0] = 0;
  translateAndExpect(xAdd, peExpect, pair);
  checkRestart(pair, "slitlj");
  addMols(pair, 300);
  for (int iAtom = 0; iAtom < space->nMol(); ++iAtom) {
    EXPECT_LT(fabs(space->x(iAtom, 0)), 5);
  }
  // pair->printXYZ("hi", 1);
}

TEST(PairFieldSlitSW, sw) {
  auto space = feasst::makeSpace({{"dimen", "2"}, {"boxLength", "12"}});
  auto pair = feasst::makePairFieldSlitSW(space);
  pair->initData("../forcefield/data.atom");
  std::vector<double> xAdd(space->dimen(), 0.);
  pair->addMol(xAdd);
  pair->initRCut(0, 1.5);
  pair->initSlit(0, 5, -5);
  pair->initEnergy();
  EXPECT_NEAR(0, pair->peTot(), feasst::DTOL);
  xAdd[0] = 3.499;
  translateAndExpect(xAdd, 0., pair);
  xAdd[0] = 3.501;
  translateAndExpect(xAdd, -1, pair);
  xAdd[0] = 4.499;
  translateAndExpect(xAdd, -1, pair);
  xAdd[0] = 4.501;
  translateAndExpect(xAdd, NUM_INF, pair);
}

TEST(PairFieldSphereSW, sw) {
  auto space = feasst::makeSpace({{"dimen", "3"}, {"boxLength", "12"}});
  auto pair = feasst::makePairFieldSphereSW(space);
  pair->initData("../forcefield/data.atom");
  std::vector<double> xAdd(space->dimen(), 0.);
  pair->addMol(xAdd);
  pair->initRCut(0, 1.5);
  pair->initSphere(5);
  pair->initEnergy();
  EXPECT_NEAR(0, pair->peTot(), feasst::DTOL);
  xAdd[0] = 3.499;
  translateAndExpect(xAdd, 0., pair);
  xAdd[0] = 3.501;
  translateAndExpect(xAdd, -1, pair);
  xAdd[0] = 4.499;
  translateAndExpect(xAdd, -1, pair);
  checkRestart(pair, "sphereSW");
  xAdd[0] = 4.501;
  translateAndExpect(xAdd, NUM_INF, pair);

  addMols(pair, 300);
  for (int iAtom = 0; iAtom < space->nMol(); ++iAtom) {
    const double dx = space->x(iAtom, 0),
                 dy = space->x(iAtom, 1),
                 dz = space->x(iAtom, 2);
    EXPECT_LT(sqrt(dx*dx + dy*dy + dz*dz), pair->radius());
  }
  // pair->printXYZ("hi", 1);
}

TEST(PairFieldSphereSW, 1dError) {
  auto space = feasst::makeSpace({{"dimen", "1"}});
  try {
    auto pair = feasst::makePairFieldSphereSW(space);
    CATCH_PHRASE("Assumes dimension > 1");
  }
}

TEST(PairFieldCylinderSW, sw) {
  auto space = feasst::makeSpace({{"dimen", "3"}, {"boxLength", "12"}});
  auto pair = feasst::makePairFieldCylinderSW(space);
  pair->initData("../forcefield/data.atom");
  std::vector<double> xAdd(space->dimen(), 0.);
  pair->addMol(xAdd);
  pair->initRCut(0, 1.5);
  const double radius = 5.;
  const int axis_dim = 2;
  pair->initCylinder(axis_dim, radius);
  pair->initEnergy();
  EXPECT_NEAR(0, pair->peTot(), feasst::DTOL);
  xAdd[0] = 3.499;
  translateAndExpect(xAdd, 0., pair);
  xAdd[0] = 3.501;
  translateAndExpect(xAdd, -1, pair);
  xAdd[0] = 4.499;
  translateAndExpect(xAdd, -1, pair);
  xAdd[axis_dim] = 4.;
  translateAndExpect(xAdd, -1, pair);
  checkRestart(pair, "cylSW");
  xAdd[1] = 5;
  translateAndExpect(xAdd, NUM_INF, pair);
  xAdd[0] = 4.501;
  xAdd[1] = 0;
  translateAndExpect(xAdd, NUM_INF, pair);

  addMols(pair, 300);
  for (int iAtom = 0; iAtom < space->nMol(); ++iAtom) {
    const double dx = space->x(iAtom, 0),
                 dy = space->x(iAtom, 1);
    EXPECT_LT(sqrt(dx*dx + dy*dy), pair->radius());
  }
  // pair->printXYZ("hi", 1);
}

TEST(PairFieldCylinderSW, 2dError) {
  auto space = feasst::makeSpace({{"dimen", "2"}});
  try {
    auto pair = feasst::makePairFieldCylinderSW(space);
    CATCH_PHRASE("assumes 3D");
  }
}

TEST(PairFieldSlitLinear, oppositeCharge) {
  const double boxLength = 8;
  const int confine_dimension = 2;
  const double upper = 2.5;
  const double lower = -2.5;
  auto space = feasst::makeSpace(
   {{"dimen", "3"},
    {"boxLength", feasst::str(boxLength)}});
  space->initBoxLength(upper*2 + 20, confine_dimension);
	auto pair1 = feasst::makePairFieldSlitLinear(space);
  pair1->initData("../forcefield/data.rpm_plus");
  pair1->initData("../forcefield/data.rpm_minus");
  pair1->initSlit(confine_dimension, upper, 2.*lower);
	pair1->initEps(0, 1.);
	pair1->initEps(1, -1.);
  pair1->initRCut(0, upper);
  pair1->initRCut(1, upper);
  feasst::PairFieldSlitLinear * pair2 = pair1->clone(space.get());
  pair2->initSlit(confine_dimension, 2*upper, lower);
	pair2->initEps(0, -1.);
	pair2->initEps(1, 1.);
	auto pair = feasst::makePairHybrid(space.get());
	pair->addPair(pair1.get());
	pair->addPair(pair2);
	pair->initEnergy();

	// test particle type 1
	std::vector<double> xAdd(space->dimen(), 0.);
  pair->addMol(xAdd, "../forcefield/data.rpm_plus");
  pair->initEnergy();
  EXPECT_NEAR(0, pair->peTot(), feasst::DTOL);
  xAdd[2] = 2.00000000001;
	translateAndExpect(xAdd, NUM_INF, pair);
  xAdd[2] = 0.5;
  translateAndExpect(xAdd, -1./5., pair);
  xAdd[2] = 2.0;
	translateAndExpect(xAdd, -4./5., pair);
  xAdd[2] = -2.0;
	translateAndExpect(xAdd, 4./5., pair);
  xAdd[2] = -2.00000000001;
	translateAndExpect(xAdd, NUM_INF, pair);

// test particle type 2
	EXPECT_EQ(1, space->nMol());
	std::vector<int> mpart = {0};
  // pair->delPart(mpart);
	space->delPart(mpart);
	EXPECT_EQ(0, space->nMol());
	xAdd[2] = 0.;
	pair->addMol(xAdd, "../forcefield/data.rpm_minus");
  pair->initEnergy();
  EXPECT_NEAR(0, pair->peTot(), feasst::DTOL);
  xAdd[2] = 2.00000000001;
	translateAndExpect(xAdd, NUM_INF, pair);
  xAdd[2] = 0.5;
  translateAndExpect(xAdd, 1./5., pair);
  xAdd[2] = 2.0;
	translateAndExpect(xAdd, 4./5., pair);
  xAdd[2] = -2.0;
	translateAndExpect(xAdd, -4./5., pair);
  xAdd[2] = -2.00000000001;
	translateAndExpect(xAdd, NUM_INF, pair);
}
