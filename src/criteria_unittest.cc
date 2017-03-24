#include <gtest/gtest.h>
#include <limits.h>
#include "criteria.h"
#include "criteria_metropolis.h"
#include "criteria_wltmmc.h"

TEST(Criteria, Metropolis) {
  const double beta = 1e-15, activ = 1;
  CriteriaMetropolis c(beta, activ);
  EXPECT_EQ(1, c.accept(0., 0., "move", 0));
}

TEST(Criteria, WLTMMC) {
  const double beta = 1e-15, activ = 1;
  const int nbin = 100;
  CriteriaWLTMMC c(beta, activ, "energy", 0, 1, nbin);
  EXPECT_EQ(0.995, c.lastbin2m());
}

TEST(Criteria, readColWLTMMC) {
  const int nMolMin = 0, nMolMax = 265;
  CriteriaWLTMMC c(0, 0, "nmol", nMolMin-0.5,nMolMax+0.5,nMolMax-nMolMin+1);
  c.readCollectMat("../unittest/colMat.txt");
  vector<int> max = findLocalMaxima(c.lnPI(), 10);
  EXPECT_EQ(2, int(max.size()));
  for (int i = 0; i < int(max.size()); ++i) {
    //cout << max[i] << " " << c.lnPI()[max[i]] << endl;
  }
}

TEST(Criteria, reweightWLTMMC) {
  const int nMolMin = 0, nMolMax = 265;
  const std::string fileName("../unittest/colMat1.txt");
  const std::string fileNameOut("colMatHWHrw.txt");
  const double temp = 525, activ = exp(-8.08564), beta = 1./(temp*8.3144621/1000);
  CriteriaWLTMMC c(beta, activ, "nmol", nMolMin-0.5,nMolMax+0.5,nMolMax-nMolMin+1);
  c.readCollectMat(fileName.c_str());
  c.findSat();
}

TEST(Criteria, reweightWLTMMCequltPatchANDnMolResizeWindow) {
  const int nMolMin = 0, nMolMax = 255;
  const std::string fileName("../unittest/colMat3.txt");
  const std::string fileNameOut("colMatasdfrw.txt");
  const double temp = 0.7, activ = 0.0183156, beta = 1./temp;
  CriteriaWLTMMC c(beta, activ, "nmol", nMolMin-0.5,nMolMax+0.5,nMolMax-nMolMin+1);
  c.readCollectMat(fileName.c_str());
  c.findSat();
  //c.lnPIrwsat(exp(-4.0));
  //EXPECT_NEAR(-8.14, log(c.activrw()), 1e-18);
  vector<int> max = findLocalMaxima(c.lnPIrw(), 10);
  c.printRWinit();
  //c.printCollectMat(fileNameOut.c_str());
  EXPECT_EQ(240, c.nMolResizeWindow(-20, 5));
  //c.printRWinit();
  //c.printCollectMat(fileNameOut.c_str());
}

TEST(Criteria, pressureANDavMacro) {
  {
    const int nMolMin = 0, nMolMax = 370;
    CriteriaWLTMMC c(1./1.5, exp(-1.568214), "nmol", nMolMin-0.5,nMolMax+0.5,nMolMax-nMolMin+1);
    c.readCollectMat("../unittest/lj/srsw/eostmmc/1.5/lj.msdb.t150.1.p_macro.dat");
    EXPECT_NEAR(310.432325, c.lnPIaverage(), 1e-6);
    EXPECT_NEAR(1., c.lnPIarea(), 1e-14);
    EXPECT_EQ(1, c.lnPInumPhases());
    vector<int> min = c.lnPIphaseBoundary();
    EXPECT_EQ(0, int(min.size()));
    //c.lnPIpressureIso(512);
    //c.printCollectMat("h123");
    vector<CriteriaWLTMMC> cvec = c.phaseSplit(c.lnPI());
    EXPECT_EQ(1, int(cvec.size()));
  }

  // lower temperature two-phase
  const int nMolMin = 0, nMolMax = 475;
  CriteriaWLTMMC c(1./0.7, exp(-5.943376), "nmol", nMolMin-0.5,nMolMax+0.5,nMolMax-nMolMin+1);
  c.readCollectMat("../unittest/lj/srsw/eostmmc/0.7/lj.msdb.t070.1.p_macro.dat");
  vector<CriteriaWLTMMC> cvec = c.phaseSplit(c.lnPI());
  EXPECT_EQ(2, c.lnPInumPhases());
  EXPECT_EQ(2, int(cvec.size()));
  EXPECT_NEAR(0., cvec[0].lnPIarea(), 1e-14);
  EXPECT_NEAR(1., cvec[1].lnPIarea(), 1e-14);
  EXPECT_NEAR(1.4205549976272893, cvec[0].lnPIaverage(), 1e-9);
  EXPECT_NEAR(4.3782625203E+02, cvec[1].lnPIaverage(), 1e-9);
  vector<int> min = c.lnPIphaseBoundary();
  EXPECT_EQ(1, int(min.size()));
  EXPECT_EQ(150, min[0]);
  EXPECT_NEAR(437.8262520297539, c.lnPIaverage(), 1e-18);
//  EXPECT_NEAR(cvec[0].lnPIpressure(512), cvec[1].lnPIpressure(512), 1e-18);
  EXPECT_NEAR(
    c.lnPIpressure(512),
    //(-c.lnPI().front() + log(cvec[0].lnPIarea()))/double(512)*0.7,
    (-c.lnPI().front() + log(cvec[1].lnPIarea()))/double(512)*0.7,
    1e-15);

  c.findSat();
  vector<CriteriaWLTMMC> cvec2 = c.phaseSplit(c.lnPIrw());
  EXPECT_EQ(215, cvec2[0].lastbin2m());
  EXPECT_EQ(216, cvec2[1].bin2m(0));
  EXPECT_NEAR(0.5, cvec2[0].lnPIarea(), 1e-6);
  EXPECT_NEAR(0.5, cvec2[1].lnPIarea(), 1e-6);
//  cout << "av " << c.lnPIrwaverage() << endl;
//  EXPECT_NEAR(
//    (-c.lnPIrw().front() + log(0.5))/double(512)/0.7,
//    0,
//    //(-c.lnPIrw().front() + log(0.5))/double(512)/0.7,
//    1e-18);
//  //EXPECT_NEAR(cvec2[0].lnPIpressure(512), cvec2[1].lnPIpressure(512), 1e-18);
  EXPECT_NEAR(1.0226696103752704, cvec2[0].lnPIaverage(), 1e-9);
  EXPECT_NEAR(431.86078547025141, cvec2[1].lnPIaverage(), 1e-9);
  c.printRWinit();
  //c.lnPIpressureIso(512);
  //c.printCollectMat("t123");
}

