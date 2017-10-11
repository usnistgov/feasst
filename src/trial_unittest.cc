/**
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of this agreement (see LICENSE.txt) and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#include <gtest/gtest.h>
#include "space.h"
#include "pair.h"
#include "pair_ideal.h"
#include "pair_lj.h"
#include "pair_lj_coul_ewald.h"
#include "functions.h"
#include "criteria.h"
#include "criteria_metropolis.h"
#include "criteria_wltmmc.h"
#include "trial.h"
#include "trial_transform.h"
#include "trial_add.h"
#include "trial_delete.h"

using namespace feasst;

int verbose_ = 0;
void vout_(std::ostream& message) { if (verbose_ == 1) { std::string data = dynamic_cast<std::ostringstream&>(message).str(); std::cout << data; message.clear(); } };

TEST(Trial, cloneANDreconstruct) {
  const double rAbove = 3, rBelow = 1;
  Space s(3, 0);
  s.init_config(12);
  s.addMolInit("../forcefield/data.atom");
  PairLJ p(&s, 5);
  p.initNeighList(rAbove, rBelow);
  p.buildNeighList();
  CriteriaMetropolis c(0.5, 0.01);
  TrialTransform tt(&p, &c, "translate");
  tt.maxMoveParam = 5;
  TrialDelete td(&p, &c);
  TrialAdd ta(&p, &c, "../forcefield/data.atom");

  ranInitByDate();
  p.initEnergy();
  const int nAttempts = 300;
  for (int i = 0; i < nAttempts; ++i) {
    td.attempt();
    tt.attempt();
  }
  const double petot = p.peTot();
  const double peLJ = p.peLJ();
  const double peLRC = p.peLRC();

  // clone, run trial, and make sure they are different
  shared_ptr<Space> s2 = s.cloneShrPtr();
  Pair* p2 = p.clone(s2.get());
  CriteriaMetropolis* c2 = c.clone();
  shared_ptr<Trial> tt2 = tt.cloneShrPtr(p2, c2);
  shared_ptr<Trial> ta2 = ta.cloneShrPtr(p2, c2);
  shared_ptr<Trial> td2 = td.cloneShrPtr(p2, c2);

  // simulate
  for (int i = 0; i < nAttempts; ++i) {
    td2->attempt();
    tt2->attempt();
    ta2->attempt();
  }

  // check that original is unchanged, and different from the duplicate
  EXPECT_NEAR(petot, p.peTot(), 1e-11);
  EXPECT_NEAR(peLJ, p.peLJ(), 1e-12);
  EXPECT_NEAR(peLRC, p.peLRC(), 1e-12);
  EXPECT_NE(petot, p2->peTot());
  EXPECT_EQ(1, p.checkNeigh());
  EXPECT_EQ(1, p2->checkNeigh());

  // ensure that random numbers in classes are different
  for (int i = 0; i < 2; ++i) {
    EXPECT_NE(s.uniformRanNum(), p.uniformRanNum());
    EXPECT_NE(p.uniformRanNum(), p.uniformRanNum());
    EXPECT_NE(c.uniformRanNum(), p.uniformRanNum());
    EXPECT_NE(tt.uniformRanNum(), p.uniformRanNum());
    EXPECT_NE(ta.uniformRanNum(), p.uniformRanNum());
    EXPECT_NE(td.uniformRanNum(), p.uniformRanNum());
  }

  // free memory
  delete p2;
  delete c2;

  // restart
  td.writeRestart("tmp/tdrst");
  TrialDelete td3("tmp/tdrst", &p, &c);
  ta.writeRestart("tmp/tarst");
  TrialAdd ta3("tmp/tarst", &p, &c);
  tt.writeRestart("tmp/ttrst");
  TrialTransform tt3("tmp/ttrst", &p, &c);
}

TEST(Trial, allmoves) {
  ranInitByDate();
  // ranInitForRepro(1491827115);
  // ranInitForRepro(1494504908);
  // declare pair and space objects for each simulation type
  Space sSPCE(3,0);
  double boxl = 24.8586887;
  for (int dim=0; dim < sSPCE.dimen(); ++dim) sSPCE.lset(boxl,dim);
  sSPCE.readXYZBulk(3, "water", "../unittest/spce/test.xyz");
  //sSPCE.readXYZBulk(3, "water", "../test/spce/test52.xyz");
  sSPCE.addMolInit("../forcefield/data.spce");
  PairLJCoulEwald pSPCE(&sSPCE, 8);
  pSPCE.initBulkSPCE(5.6, 38);

  Space sLJ(3,1);
  sLJ.init_config(12);
  sLJ.addMolInit("../forcefield/data.atom");
  PairLJ pLJ(&sLJ, 5);

  Space sID(3,2);
  sID.init_config(12);
  sID.addMolInit("../forcefield/data.atom");
  PairIdeal pID(&sID, 5);

  // loop through each pair type listed in string vector pairType
  vector<std::string> pairType;
  pairType.push_back("id");
  pairType.push_back("lj");
  //pairType.push_back("spce");
  for (vector<std::string>::iterator pt = pairType.begin(); pt != pairType.end(); ++pt) {

    double rAbove=0, rBelow=0, beta=0, activ=0, maxMoveParam=0;
    int nAttempts = 0;
    std::string addType;
    Pair* p = NULL;
    Space* s = NULL;
    if ((*pt).compare("lj") == 0) {
      p = &pLJ;
      s = &sLJ;
      nAttempts = 50;
      rAbove = 3, rBelow = 1, beta = 1, activ = 0.01, maxMoveParam = 5;
      addType.assign("../forcefield/data.atom");
    } else if ((*pt).compare("spce") == 0) {
      p = &pSPCE;
      s = &sSPCE;
      nAttempts = 5;
      rAbove = 5., rBelow = 2.5, beta = 0.25, activ = 0.01, maxMoveParam = 0.1;
      addType.assign("../forcefield/data.spce");
    } else if ((*pt).compare("id") == 0) {
      p = &pID;
      s = &sID;
      nAttempts = 50;
      rAbove = 3, rBelow = 1, beta = 0.025, activ = 0.01, maxMoveParam = sID.minl()/2.;
      addType.assign("../forcefield/data.atom");
    }

    // loop through criteria types
    vector<std::string> critType;
    critType.push_back("metropolis");
    critType.push_back("wltmmc");
    for (vector<std::string>::iterator ct = critType.begin(); ct != critType.end(); ++ct) {
      std::cout << *pt << " " << *ct << std::endl;
      CriteriaMetropolis cm(beta, activ);
      CriteriaWLTMMC cwltmmc(beta, activ, "nmol", 0 - 0.5, (*s).nMol() + 0.5, (*s).nMol() + 1);
      //CriteriaWLTMMC c(beta, activ, "energy", -300, -250, 100);
      Criteria *c = NULL;
      if ((*ct).compare("metropolis") == 0) {
        c = &cm;
      } else if ((*ct).compare("wltmmc") == 0) {
        c = &cwltmmc;
      }
      (*p).initNeighList(rAbove, rBelow);
      (*p).initEnergy();
      TrialTransform tt(p,c, "translate");
      tt.maxMoveParam = maxMoveParam;
      TrialTransform tr(p,c, "rotate");
      TrialTransform tst(p,c, "smctrans");
      TrialAdd ta(p,c, addType.c_str());
      TrialAdd tamfb(p,c, addType.c_str());
      tamfb.numFirstBeads(3);
      TrialAdd taavb(p,c, addType.c_str());
      taavb.initAVB(rAbove, rBelow);
      TrialAdd taavbmfb(p,c, addType.c_str());
      taavbmfb.numFirstBeads(3);
      taavbmfb.initAVB(rAbove, rBelow);
      TrialDelete td(p,c);
      TrialDelete tdmfb(p,c);
      tdmfb.numFirstBeads(3);
      TrialDelete tdavb(p,c);
      tdavb.initAVB(rAbove, rBelow);
      TrialDelete tdavbmfb(p,c);
      tdavbmfb.numFirstBeads(3);
      tdavbmfb.initAVB(rAbove, rBelow);

      for (int i = 0; i < nAttempts; ++i) {
        vout_(std::ostringstream().flush() << "attempting del " << i << " " << (*s).nMol() << std::endl );
        td.attempt();
        vout_(std::ostringstream().flush() << "attempting trans " << i << " " << (*s).nMol() << std::endl );
        tt.attempt();
        vout_(std::ostringstream().flush() << "attempting rotate " << i << " " << (*s).nMol() << std::endl );
        tr.attempt();
        vout_(std::ostringstream().flush() << "attempting insert " << i << " " << (*s).nMol() << std::endl );
        ta.attempt();
        vout_(std::ostringstream().flush() << "attempting deletion mfb" << i << " " << (*s).nMol() << std::endl );
        tdmfb.attempt();
        vout_(std::ostringstream().flush() << "attempting avb deletion mfb" << i << " " << (*s).nMol() << std::endl );
        tdavbmfb.attempt();
        vout_(std::ostringstream().flush() << "attempting insert mfb" << i << " " << (*s).nMol() << std::endl );
        tamfb.attempt();
        vout_(std::ostringstream().flush() << "attempting avb insert mfb" << i << " " << (*s).nMol() << std::endl );
        taavbmfb.attempt();
        vout_(std::ostringstream().flush() << "attempting avb " << i << " " << (*s).nMol() << std::endl );
        taavb.attempt();
        vout_(std::ostringstream().flush() << "attempting avb deletion " << i << " " << (*s).nMol() << std::endl );
        tdavb.attempt();
        vout_(std::ostringstream().flush() << "attempting avbmfb " << i << " " << (*s).nMol() << std::endl );
      }

      const double petot = (*p).peTot();
      double peLJ=0, peLRC=0, peQReal=0, peQFrr=0, peQFrrSelf=0;
      if ((*pt).compare("id") == 0) {
        EXPECT_NEAR(0, petot, 1e-15);
      } else if ((*pt).compare("lj") == 0) {
        peLJ = pLJ.peLJ();
        peLRC = pLJ.peLRC();
      } else if ((*pt).compare("spce") == 0) {
        peLJ = pSPCE.peLJ();
        peLRC = pSPCE.peLRC();
        peQReal = pSPCE.peQReal();
        peQFrr = pSPCE.peQFrr();
        peQFrrSelf = pSPCE.peQFrrSelf();
      }
      (*p).initEnergy();
      EXPECT_NEAR(petot, (*p).peTot(), 1e-9);
      if ((*pt).compare("lj") == 0) {
        EXPECT_NEAR(peLJ, pLJ.peLJ(), 1e-12);
        EXPECT_NEAR(peLRC, pLJ.peLRC(), 1e-12);
      } else if ((*pt).compare("spce") == 0) {
        EXPECT_NEAR(peLJ, pSPCE.peLJ(), 1e-9);
        EXPECT_NEAR(peLRC, pSPCE.peLRC(), 1e-9);
        EXPECT_NEAR(peQReal, pSPCE.peQReal(), 1e-9);
        EXPECT_NEAR(peQFrr, pSPCE.peQFrr(), 1e-9);
        EXPECT_NEAR(peQFrrSelf, pSPCE.peQFrrSelf(), 1e-9);
        (*s).checkBond("spce", 1e-12);
      }
//      EXPECT_NE(td.acceptPer(), 0);
//      EXPECT_NE(ta.acceptPer(), 0);
//      EXPECT_NE(tt.acceptPer(), 0);
//      EXPECT_NE(tr.acceptPer(), 0);
      if ( ((*pt).compare("lj") == 0) && ((*ct).compare("metropolis") == 0) ) {
        EXPECT_EQ(tr.acceptPer(), 1);
      }
      EXPECT_EQ(td.attempted(), nAttempts);
      EXPECT_EQ(ta.attempted(), nAttempts);
      EXPECT_EQ(tt.attempted(), nAttempts);
      EXPECT_EQ(tr.attempted(), nAttempts);
//      if (pt->compare("spce") == 0) {
//        EXPECT_NE(tr.accepted(), nAttempts);
//      }
      EXPECT_EQ(1, (*p).checkNeigh());

      // test replace and restore
      if (ct->compare("wltmmc") == 0) {
        cwltmmc.lnPIupdate();
        long double lnpif = cwltmmc.lnPI().front();
        CriteriaMetropolis cnew(1, 1);
        tt.replaceCriteria(&cnew);
        for (int i = 0; i < 10; ++i) tt.attempt();
        cwltmmc.lnPIupdate();
        tt.restoreCriteria();
        EXPECT_EQ(lnpif, cwltmmc.lnPI().front());
        EXPECT_EQ(tt.criteria(), &cwltmmc);
      }
    }
  }
}

