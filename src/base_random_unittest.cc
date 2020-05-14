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
#include <limits.h>
#include "base_random.h"
#include "accumulator.h"
#include "histogram.h"

using namespace feasst;

TEST(BaseRandom, uniformRanNum) {
  ranInitByDate();
  BaseRandom ran;
  const int min = 5, max = 8, n = 10000;
  vector<int> x(max);
  for (int i = 0; i < n; ++i) {
    ++x.at(ran.uniformRanNum(min, max) - 1);
  }
  for (int i = min; i <= max; ++i) {
    EXPECT_NEAR(1./(max - min + 1), double(x.at(i-1))/n, 5e-2);
  }
  EXPECT_EQ(0, x.at(min-2));

//  // check if the restart is the same random number sequence.
//  ran.writeRngRestart("tmp/rstnr3");
//  BaseRandom ran2("tmp/rstnr3");
//  EXPECT_NEAR(ran.uniformRanNum(), ran2.uniformRanNum(), 1e-19);
}

TEST(BaseRandom, stdNormRanNum) {
  ranInitByDate();
  BaseRandom ran;
  const int n = 10000;
  Accumulator a;
  for (int i = 0; i < n; ++i) {
    const double random = ran.stdNormRanNum();
    a.accumulate(random);
  }
  EXPECT_NEAR(0, a.average(), 5e-2);
  EXPECT_NEAR(1, a.std(), 5e-2);

}

TEST(BaseRandom, gaussRanNum) {
  ranInitByDate();
  BaseRandom ran;
  const int n = 10000;
  const double sig = 5, av = 10;
  Accumulator a;
  Histogram hist(0.5);
  for (int i = 0; i < n; ++i) {
    const double random = ran.gaussRanNum(sig, av);
    a.accumulate(random);
    hist.accumulate(random);
  }
  EXPECT_NEAR(av, a.average(), sig*5e-2);
  EXPECT_NEAR(sig, a.std(), sig*5e-2);
  hist.write("tmp/histhist");
}

TEST(BaseRandom, hash) {
  BaseRandom ran;
  std::string hash = ran.randomHash();
  EXPECT_NE(hash, ran.randomHash());
  EXPECT_NE(hash, ran.randomHash());
}

TEST(BaseRandom, ranShell) {
  BaseRandom math;
  ranInitByDate();
  const double rabove = 1.08, rbelow = 1.0;

  for (int dim = 2; dim < 4; ++dim) {
    for (int i = 0; i < 1000; ++i) {
      vector<double> x = math.ranShell(rabove, rbelow, dim);
//      if (dim == 2) cout << vec2str(x) << endl;
      const double r = sqrt(vecDotProd(x, x));
      EXPECT_GE(r, rbelow);
      EXPECT_LE(r, rabove);
    }
  }

  try {
    math.ranShell(1., 0.5, 1);
    CATCH_PHRASE("ranShell not implemented for current dimen");
  }
  try {
    math.ranShell(1., 0.5, 4);
    CATCH_PHRASE("ranShell not implemented for current dimen");
  }
}

TEST(BaseRandom, ranFromCPDF) {
  BaseRandom math;
  vector<double> cpdf;
  const int ncpdf = 10, n = 100;
  for (int i = 0; i < ncpdf; ++i) cpdf.push_back((i+1)/double(ncpdf));
  vector<double> cpdfran(ncpdf);
  for (int i = 0; i < n; ++i) {
    const int j = math.ranFromCPDF(cpdf);
    ++cpdfran[j];
  }
  for (int i = 0; i < ncpdf; ++i) EXPECT_NEAR(cpdfran[i]/double(n), ncpdf/double(n), 0.2);
}

TEST(BaseRandom, quatRandomNorm) {
  ranInitByDate();
  BaseRandom math;
  const int nRand = 100;
  vector<double> q = math.quatRandom(0);
  EXPECT_NEAR(0, q[0], 1e-18);
  EXPECT_NEAR(0, q[1], 1e-18);
  EXPECT_NEAR(0, q[2], 1e-18);
  EXPECT_NEAR(1, q[3], 1e-18);
  for (int i = 0; i < nRand; ++i) {
    q = math.quatRandom(1);
    EXPECT_NEAR(1, vecDotProd(q, q), 1e-15);
    vector<double> q2 = math.quatRandom();
    EXPECT_NEAR(1, vecDotProd(q2, q2), 1e-14);
  }
}

TEST(BaseRandom, randQuatAndRotSphere) {
  BaseRandom math;
  vector<double> x;
  x.resize(3);
  vector<double> y;
  y.resize(3);
  x[0] = 1;
  x[1] = 1;
  x[2] = 1;
  vector<vector<double> > r;
  const int dim1 = 3; const int dim2 = 3;
  r.resize(dim1, vector<double>(dim2));
  ranInitByDate();
  const int nRand = 1000;
  for (int i = 0; i < nRand; ++i) {
    vector<double> q = math.quatRandom();
    EXPECT_NEAR(1, q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3], 1e-14);
    //quat2rot(q, r);
    //matVecMul(r, x, y);
  }
}

TEST(BaseRandom, ranUnitSphere) {
  ranInitByDate();
  BaseRandom math;
  for (int i = 0; i < 100; ++i) {
    vector<double> x = math.ranUnitSphere(3);
    EXPECT_NEAR(1, vecDotProd(x, x), 1e-14);
  }
}

TEST(BaseRandom, ranAngle) {
  BaseRandom b;
  EXPECT_NEAR(PI/3, b.ranAngle(100, PI/3), 0.4);
  //for (int i = 0; i < 100; ++i) cout << b.ranAngle(100, PI/3) << endl;
}

TEST(BaseRandom, eulerRandom) {
  ranInitByDate();
  BaseRandom math;
  const int n = 10000;
  vector<double> xref(3); xref[0] = 1;
  for (int i = 0; i < n; ++i) {
    vector<double> euler = math.eulerRandom();
    vector<vector<double> > R = Euler2RotMat(euler);
    vector<double> x = matVecMul(R, xref);
    // test by visual inspection
    //cout << x[0] << " " << x[1] << " " << x[2] << endl;
  }
}

