#include <gtest/gtest.h>
#include "pair_ideal.h"
#include "pair_hybrid.h"
#include "pair_wall.h"
#include "pair_hs.h"
#include "pair_patch_kf.h"
#include "pair_lj.h"
#include "pair_lj_multi.h"
#include "pair_lj_coul_ewald.h"
#include "pair_tabular_1d.h"
#include "pair_hard_circle.h"
#include "pair_squarewell.h"
#include "mc_wltmmc.h"
#include "ui_abbreviated.h"
#include "trial_add.h"
#include "trial_delete.h"
#include "trial_md.h"
#include "trial_transform.h"
#include "criteria_mayer.h"

using namespace feasst;

// check that you obtain the same tmmc probability distribution function whether or not you use special moves
TEST(MC, WLTMMC_Ideal) {
  //ranInitForRepro(1494446668);
  ranInitByDate();
  // to make acceptance equal at midpoint density, (activ/rho)^2 == 1
  Space s(3, 0);
  const int nMolMax = 4, nMolMin = 1, nAttemptsSimple = 600, nAttempts = 600, ncfreq = 1;
  // these values were used in testing cell list const int nMolMax = 4, nMolMin = 1, nAttemptsSimple = 600000, nAttempts = 600000, ncfreq = 1;
  // ditto const double boxl = 12, beta = 1, rhoMid = double(nMolMax)/2/pow(boxl,s.dimen()), activ = rhoMid, rCut = 3, rAbove = 3, rBelow = 1;
  const double boxl = 6, beta = 1, rhoMid = double(nMolMax)/2/pow(boxl,s.dimen()), activ = rhoMid, rCut = 3, rAbove = 3, rBelow = 1;
  for (int dim=0; dim < s.dimen(); ++dim) s.lset(boxl,dim);
  s.addMolInit("../forcefield/data.atom");
  for (int i = 0; i < nMolMin; ++i) s.addMol("../forcefield/data.atom");
  PairIdeal p(&s, rCut);
  CriteriaWLTMMC c(beta, activ,"nmol",nMolMin-0.5,nMolMax+0.5,nMolMax-nMolMin+1);
  c.collectInit();
  c.tmmcInit();
  MC mc(&s,&p,&c);
  c.prefilColMat(1e-15);

  transformTrial(&mc, "translate");
  deleteTrial(&mc);
  addTrial(&mc, "../forcefield/data.atom");

  for (int i = 0; i < nAttemptsSimple; ++i) {
    mc.attemptTrial();
    if ( ((i+1)%ncfreq == 0) || (i+1 == nAttemptsSimple) ) {
      c.lnPIupdate();
      c.printCollectMat("tmp/coltmps.txt");
    }
  }
  shared_ptr<CriteriaWLTMMC> cSimple = c.cloneShrPtr();

  for (int t = 0; t < 3; ++t) {
    cout << "t " << t << std::endl;
    c.zeroStat();
    c.collectInit();
    c.tmmcInit();
    MC mc2(&s,&p,&c);
    c.prefilColMat(1e-15);
    if (t == 0) {
      transformTrial(&mc, "translate");

      TrialDelete td;
      td.numFirstBeads(3);
      mc2.initTrial(&td);

      TrialAdd ta("../forcefield/data.atom");
      ta.numFirstBeads(3);
      mc2.initTrial(&ta);

      EXPECT_FALSE(p.neighOn());
    } else if (t == 1) {


      transformTrial(&mc2, "translate");

      TrialDelete td(&s, &p, &c);
      td.initAVB(rAbove, rBelow);
      mc2.initTrial(&td);
      mc2.neighAVBInit(rAbove, rBelow);

      //shared_ptr<TrialDelete> td = make_shared<TrialDelete>();
      //mc2.initTrial(td);
      //td->initAVB(rAbove, rBelow);
      //mc2.neighAVBInit(rAbove, rBelow);

      //mc2.addAVBTrial(rAbove, rBelow, "../forcefield/data.atom");
      shared_ptr<TrialAdd> ta = make_shared<TrialAdd>("../forcefield/data.atom");
      mc2.initTrial(ta);
      ta->initAVB(rAbove, rBelow);
      mc2.neighAVBInit(rAbove, rBelow);

      EXPECT_TRUE(p.neighOn());
      EXPECT_EQ(0, s.cellType());
    } else if (t == 2) {
      s.updateCells(rCut, rCut);

      transformTrial(&mc2, "translate");

      //mc2.deleteAVBMFBTrial(rAbove, rBelow, 3);
      shared_ptr<TrialDelete> td = make_shared<TrialDelete>();
      mc2.initTrial(td);
      td->initAVB(rAbove, rBelow);
      td->numFirstBeads(3);
      mc2.neighAVBInit(rAbove, rBelow);

      //mc2.addAVBMFBTrial(rAbove, rBelow, "../forcefield/data.atom", 3);
      shared_ptr<TrialAdd> ta = make_shared<TrialAdd>("../forcefield/data.atom");
      mc2.initTrial(ta);
      ta->initAVB(rAbove, rBelow);
      ta->numFirstBeads(3);
      mc2.neighAVBInit(rAbove, rBelow);

      EXPECT_TRUE(p.neighOn());
      // ditto EXPECT_EQ(1, s.cellType());
    }

    for (int i = 0; i < nAttempts; ++i) {
      mc2.attemptTrial();
      if ( ((i+1)%ncfreq == 0) || (i+1 == nAttempts) ) {
        c.lnPIupdate();
        if (t == 0) c.printCollectMat("tmp/coltmp0.txt");
        if (t == 1) c.printCollectMat("tmp/coltmp1.txt");
        if (t == 2) c.printCollectMat("tmp/coltmp2.txt");
        if (t == 3) c.printCollectMat("tmp/coltmp3.txt");
      }
    }

    // expect that lnPI is within % of the simple method
    for (int i = 0; i < int(c.lnPI().size()); ++i) {
      EXPECT_NEAR(c.lnPI()[i], cSimple->lnPI()[i], 3.);
      //EXPECT_NEAR(1, c.lnPI()[i]/cSimple->lnPI()[i], 0.9);
    }

    // expect that each trial is attempted the expected number of times, with 25% statistical error
    for (int i = 0; i < mc2.nTrials(); ++i) {
      EXPECT_LT(nAttempts/mc2.nTrials()/1.25, mc2.trialVec()[i]->attempted());
    }
    if (p.neighOn()) EXPECT_EQ(1, p.checkNeigh());
  }
  EXPECT_EQ(1, p.checkEnergy(1e-11, 0));
}

TEST(MC, ljmuvtmetropANDclone) {
  const double beta = 1./2., activ = 0.97747, rCut = 2.5, boxl = pow(250., 1./3.);
  ranInitByDate();
  Space s(3,0);
  for (int dim=0; dim < s.dimen(); ++dim) s.lset(boxl,dim);
  s.addMolInit("../forcefield/data.atom");
  PairLJ p(&s, rCut);
  p.initEnergy();
  CriteriaMetropolis c(beta, activ);
  MC mc(&s,&p,&c);

  transformTrial(&mc, "translate");
  mc.weight=0.5;
  deleteTrial(&mc);
  addTrial(&mc, "../forcefield/data.atom");
  //mc.initTrial(new TrialAdd("../forcefield/data.atom"));
  //mc.addTrial("../forcefield/data.atom");
//  TrialDelete td;
//  mc.initTrial(&td);
//
//  TrialAdd ta("../forcefield/data.atom");
//  mc.initTrial(&ta);

  // attempt some monte carlo trials
  mc.initLog("tmp/ljloglog", 1e3);
  const int nAttempts = 3000;
  for (int i = 0; i < nAttempts; ++i) {
    mc.attemptTrial();
  }

  // expect that each trial is attempted the expected number of times, with 15% statistical error
  EXPECT_LT(nAttempts/2/1.15, mc.trialVec()[0]->attempted());
  EXPECT_LT(nAttempts/4/1.15, mc.trialVec()[1]->attempted());
  EXPECT_LT(nAttempts/4/1.15, mc.trialVec()[2]->attempted());

  EXPECT_EQ(1, p.checkEnergy(1e-11, 0));

  // clone simulation, run trials in clone, and expect no change in the original
  double petot = mc.pePerMol();
  shared_ptr<MC> mc2 = mc.cloneShrPtr();
  for (int i = 0; i < nAttempts; ++i) mc2->attemptTrial();
  EXPECT_EQ(petot, mc.pePerMol());
  if (petot*mc2->pePerMol()!=0) EXPECT_NE(petot, mc2->pePerMol());

  // restart simulation from file, run trials, and expect no change in the original
  petot = mc.pePerMol();
  mc.writeRestart("tmp/ljrst");
  MC mc3("tmp/ljrst");
  for (int i = 0; i < nAttempts; ++i) mc3.attemptTrial();
  EXPECT_EQ(petot, mc.pePerMol());
  if (petot*mc3.pePerMol()!=0) {
    EXPECT_NE(petot, mc3.pePerMol());
  }
}

TEST(MC, ljnvtmetropANDremoveTrial) {
  const double beta = 1/2, rCut = 2.5;
  ranInitByDate();
  Space s(3,0);
  s.init_config(12);
  PairLJ p(&s, rCut);
  p.buildNeighList();
  p.initEnergy();
  CriteriaMetropolis c(beta, 0);
  MC mc(&s,&p,&c);

  transformTrial(&mc, "translate");

  // attempt some monte carlo trials
  const int nAttempts = 300;
  for (int i = 0; i < nAttempts; ++i) {
    mc.attemptTrial();
    //std::cout << p.peTot() << std::endl;
  }

  EXPECT_EQ(1, p.checkEnergy(1e-7, 0));
  
  transformTrial(&mc, "gca");
  mc.removeTrial(0);
  EXPECT_EQ(1, mc.nTrials());
  for (int i = 0; i < nAttempts; ++i) {
    mc.attemptTrial();
  }
}

TEST(MC, ljmuvttmmc) {
  const double rCut = 3., beta = 1./1.5, activ = exp(-1.568214), boxl = pow(512, 1./3.);
  const int nMolMax = 3;
  ranInitByDate();
  Space s(3,0);
  for (int dim=0; dim < s.dimen(); ++dim) s.lset(boxl,dim);
  s.addMolInit("../forcefield/data.atom");
  PairLJ p(&s, rCut);
  p.initEnergy();
  CriteriaWLTMMC c(beta, activ,"nmol0",0-0.5,nMolMax+0.5,nMolMax+1);
  MC mc(&s,&p,&c);

  transformTrial(&mc, "translate");
  mc.weight=0.5;
  deleteTrial(&mc);
  addTrial(&mc, "../forcefield/data.atom");

  // attempt some monte carlo trials
//  for (int i = 0; i < 3000; ++i) {
//    mc.attemptTrial();
//    std::cout << i << " " << s.nMol() << " " << p.peTot()/s.nMol() << std::endl;
//  }
  c.collectInit();
//  for (int i = 0; i < 3000; ++i) {
//    mc.attemptTrial();
//    std::cout << i << " " << s.nMol() << " " << p.peTot()/s.nMol() << std::endl;
//  }
  c.tmmcInit();
  for (int i = 0; i < 500; ++i) {
    mc.attemptTrial();
//    cout << i << " " << s.nMol() << " " << mc.pePerMol() << " " << c.nSweep() << endl;
//    { std::ostringstream s;
//      s << "tmp" << i << ".txt";
//      c.printCollectMat(s.str().c_str());
//    }
  }
  c.lnPIupdate();
  c.printCollectMat("tmp/tmp.txt");
  EXPECT_EQ(1, p.checkEnergy(1e-9, 0));
}

TEST(MC, muvttmmcspce) {
  const double boxl = 20., rCut = boxl/2., temp = 525, activ = exp(-8.08564), beta = 1./(temp*8.3144621/1000);
  const int nMolMin = 0, nMolMax = 20, ncfreq = 100, npr = 250, nprPerMacrostate = npr/(nMolMax - nMolMin + 1) + 1;

  // initialize
  ranInitByDate();
  const int nThreads = 2;
  vector<vector<int> > nMolVec = nWindow(0, nMolMax, 1.33, nThreads, 2);
  vector<int>nprv;
  for (int i = 0; i < int(nMolVec.size()); ++i) {
    nprv.push_back((nMolVec[i][1] - nMolVec[i][0] + 1) * nprPerMacrostate);
    //std::cout << "# MC " << i << " " << nprv.back() << std::endl;
  }
  vector<shared_ptr<Space> > s(nThreads);
  vector<shared_ptr<PairLJCoulEwald> > p(nThreads);
  vector<shared_ptr<CriteriaWLTMMC> > c(nThreads);
  vector<shared_ptr<WLTMMC> > mc(nThreads);
  for (int t = 0; t < nThreads; ++t) {
    s[t] = make_shared<Space>(3,t);
    for (int dim=0; dim < s[t]->dimen(); ++dim) s[t]->lset(boxl,dim);
    s[t]->addMolInit("../forcefield/data.spce");   // add one molecule in order to initialize ntype array
    p[t] = make_shared<PairLJCoulEwald>(s[t].get(), rCut);
    p[t]->initBulkSPCE(5.6, 38);
    c[t] = make_shared<CriteriaWLTMMC>(beta, activ,"nmol",nMolVec[t][0]-0.5,nMolVec[t][1]+0.5,nMolVec[t][1]-nMolVec[t][0]+1);
    mc[t] = make_shared<WLTMMC>(s[t].get(), p[t].get(), c[t].get());

    transformTrial(mc[t], "translate");
    transformTrial(mc[t], "rotate");
    deleteTrial(mc[t]);
    addTrial(mc[t], "../forcefield/data.spce");

    TrialDelete tdmfb;
    tdmfb.numFirstBeads(10);
    mc[t]->initTrial(&tdmfb);

    TrialAdd tamfb("../forcefield/data.spce");
    tamfb.numFirstBeads(10);
    mc[t]->initTrial(&tamfb);
    mc[t]->nMolSeekInRange();
    c[t]->collectInit();
    c[t]->tmmcInit();
    for (long long i = 0; i < 250; ++i) {
      //std::cout << i << " " << s.nMol() << std::endl;
      mc[t]->attemptTrial();
      if (i % ncfreq == 0) {
        c[t]->lnPIupdate();
        EXPECT_EQ(1, p[t]->checkEnergy(1e-11, 0));
      }
    }
  }
//  c.front()->printCollectMat("tmp/colMat.txt", c);
}

TEST(MC, nseek) {
  const double boxl = 20., rCut = boxl/2., temp = 525, activ = exp(-8.08564), beta = 1./(temp*8.3144621/1000);
  const int nMolMax = 20;

  // initialize
  Space s(3,0);
  for (int dim=0; dim < s.dimen(); ++dim) s.lset(boxl,dim);
  s.addMolInit("../forcefield/data.spce");   // add one molecule in order to initialize ntype array
  PairLJCoulEwald p(&s, rCut);
  p.initBulkSPCE(5.6, 38);
  CriteriaWLTMMC c(beta, activ,"nmol",0-0.5,nMolMax+0.5,nMolMax+1);
  MC mc(&s,&p,&c);
  transformTrial(&mc, "translate");
  transformTrial(&mc, "rotate");
  mc.nMolSeek(20, "../forcefield/data.spce", 1e5);
  EXPECT_EQ(20, s.nMol());
  mc.nMolSeek(2, "../forcefield/data.spce", 1e5);
  EXPECT_EQ(2, s.nMol());
}

TEST(MC, equltl43muvttmmcANDinitWindows) {
  const double temp = 1., activ = exp(-2.), boxl = 9, beta = 1/temp;
  const int nMolMax = 50, nMolMin = 10, npr = 200;
  Space s(3, 0);
  for (int dim=0; dim < s.dimen(); ++dim) s.lset(boxl,dim);
  s.addMolInit("../forcefield/data.equltl43");
  PairLJ p(&s, pow(2, 1./6.));
  p.initLMPData("../forcefield/data.equltl43");
  p.cutShift(1);
  p.lrcFlag = 0;
  p.checkEnergy(1e-9, 1);
//  s.updateCells(p.rCut(), p.rCut());
  CriteriaWLTMMC c(beta, activ,"nmol",nMolMin-0.5,nMolMax+0.5,nMolMax-nMolMin+1);
  WLTMMC mc(&s,&p,&c);

  transformTrial(&mc, "translate");
  transformTrial(&mc, "rotate");
  deleteTrial(&mc);
  addTrial(&mc, "../forcefield/data.equltl43");
 // mc.initTrial(new TrialDelete());
 // mc.initTrial(new TrialAdd("../forcefield/data.equltl43"));
//  mc.deleteTrial();
//  mc.addTrial("../forcefield/data.equltl43");

//  TrialTransform tt("translate");
//  mc.initTrial(&tt);
//  TrialTransform tr("rotate");
//  mc.initTrial(&tr);
//
//  TrialDelete td;
//  mc.initTrial(&td);
//
//  TrialAdd ta("../forcefield/data.equltl43");
//  mc.initTrial(&ta);

  c.collectInit();
  c.tmmcInit();
  
  // if this test segfaults, check pointer ownership
  // try commenting out the windows lines to see if test works
  // its possible the number of processors is too large for 
  //  the number of particles (nMolMax)
  mc.initWindows(1);
  
  mc.setNFreqCheckE(npr/2, 1e-9);
  mc.initColMat("tmp/coll", 2*npr);
  mc.initLog("tmp/ll", 1e3);
  mc.initRestart("tmp/llr", 2*npr);
  mc.runNumSweeps(0, npr);
//  mc.nMolSeekInRange();
//  for (long long i = 0; i < npr; ++i) {
//    mc.attemptTrial();
//    //cout << s.nMol() << " " << mc.pePerMol() << endl;
//    p.checkEnergy(1e-9, 0);
//    if (i%100 == 0) {
//      c.lnPIupdate();
//    }
//  }
//  p.checkEnergy(1e-9, 0);
}

TEST(MC, wltmmccloneANDreconstruct) {
  const int nMolMax = 100, npr = 200;
  Space s(3, 0);
  for (int dim=0; dim < s.dimen(); ++dim) s.lset(9,dim);
  s.readXYZBulk(4, "../forcefield/data.equltl43", "../unittest/equltl43/two.xyz");
  s.addMolInit("../forcefield/data.equltl43");
  PairLJ p(&s, pow(2, 1./6.));
  p.initLMPData("../forcefield/data.equltl43");
  p.cutShift(1);
  p.lrcFlag = 0;
  p.checkEnergy(1e-9, 1);
  CriteriaWLTMMC c(1., exp(-2.), "nmol",-0.5,nMolMax+0.5,nMolMax+1);
  WLTMMC mc(&s,&p,&c);

  transformTrial(&mc, "translate");
  transformTrial(&mc, "rotate");
  deleteTrial(&mc);
  addTrial(&mc, "../forcefield/data.equltl43");
  //mc.initTrial(new TrialDelete());
  //mc.initTrial(new TrialAdd("../forcefield/data.equltl43"));

  //mc.transformTrial("translate");
  //mc.transformTrial("rotate");
  //mc.deleteTrial();
  //mc.addTrial("../forcefield/data.equltl43");
  c.collectInit();
  c.tmmcInit();
  for (long long i = 0; i < npr; ++i) {
    mc.attemptTrial();
    //cout << s.nMol() << " " << mc.pePerMol() << endl;
    p.checkEnergy(1e-9, 0);
//    cout << "peto " << p.peTot() << endl;
    if (i%100 == 0) {
      c.lnPIupdate();
    }
  }
  p.checkEnergy(1e-9, 0);
//  cout << " here 1" << endl;
  // clone simulation, run trials in clone, and expect no change in the original
  const double petot = mc.pePerMol();
//  cout << " here 2" << endl;
  shared_ptr<WLTMMC> mc2 = mc.cloneShrPtr();
//  cout << " here 2-3" << endl;
  for (int i = 0; i < npr; ++i) mc2->attemptTrial();
//  cout << " here 3" << endl;
  EXPECT_EQ(petot, mc.pePerMol());
//  cout << " here 4" << endl;
  if (petot*mc2->pePerMol() != 0) EXPECT_NE(petot, mc2->pePerMol());
//  cout << " here 5" << endl;

  EXPECT_EQ(1, mc.checkTrialCriteria());
}

//TEST(MC, mchybrid) {
//  Space s(3, 0);
//  for (int dim=0; dim < s.dimen(); ++dim) s.lset(8,dim);
//  s.addMolInit("../forcefield/data.cg3_60_43_1");
//  s.updateCells(3.);
//  PairLJMulti pLJ(&s, 3);
//  pLJ.initLMPData("../forcefield/data.cg3_60_43_1");
//  pLJ.epsijset(1, 2, 0.);
//  pLJ.epsijset(2, 1, 0.);
//  pLJ.epsijset(2, 2, 0.);
//  pLJ.linearShift(1);
//  pLJ.initEnergy();
//  PairLJMulti pWCA(&s, pow(2, 1./6.));
//  pWCA.initLMPData("../forcefield/data.cg3_60_43_1");
//  pWCA.epsijset(1, 1, 0.);
//  pWCA.cutShift(1);
//  pWCA.initEnergy();
//  PairHybrid p(&s, s.minl()/2.);
//  p.addPair(&pLJ);
//  p.addPair(&pWCA);
//  //p.initEnergy();
//  EXPECT_EQ(2, p.nPairs());
//  const int nMolMax = 50;
//  CriteriaWLTMMC c(1., exp(-2.), "nmol",-0.5,nMolMax+0.5,nMolMax+1);
//  WLTMMC mc(&s,&p,&c);
//  mc.transformTrial("translate");
//  mc.transformTrial("rotate");
//  mc.deleteTrial();
//  mc.addTrial("../forcefield/data.cg3_60_43_1");
//  c.collectInit();
//  c.tmmcInit();
//  mc.setNFreqCheckE(10, 1e-6);
//  mc.initRestart("tmp/rst", 10);
//  mc.nMolSeek(nMolMax, "../forcefield/data.cg3_60_43_1", 1e6);
//  mc.runNumTrials(100);
//  mc.nMolSeek(nMolMax, "../forcefield/data.cg3_60_43_1", 1e6);
//  EXPECT_EQ(1, mc.checkTrialCriteria());
//
//  // new pair style in one instead of hybrid
//  PairLJMulti pLJnew(&s, 3);
//  pLJnew.initLMPData("../forcefield/data.cg3_60_43_1");
//  pLJnew.rCutijset(0, 0, 0.);
//  pLJnew.rCutijset(0, 1, 0.);
//  pLJnew.rCutijset(0, 2, 0.);
//  pLJnew.rCutijset(1, 1, 3.);
//  pLJnew.linearShiftijset(1, 1, 1);
//  pLJnew.rCutijset(1, 2, pow(2, 1./6.));
//  pLJnew.cutShiftijset(1, 2, 1);
//  pLJnew.rCutijset(2, 2, pow(2, 1./6.));
//  pLJnew.cutShiftijset(2, 2, 1);
//  pLJnew.initEnergy();
//  EXPECT_NEAR(pLJnew.peTot(), p.peTot(), 1e-11);
//  CriteriaWLTMMC cnew(1., exp(-2.), "nmol",-0.5,nMolMax+0.5,nMolMax+1);
//  WLTMMC mcnew(&s,&pLJnew,&cnew);
//  mcnew.transformTrial("translate");
//  mcnew.transformTrial("rotate");
//  mcnew.deleteTrial();
//  mcnew.addTrial("../forcefield/data.cg3_60_43_1");
//  cnew.collectInit();
//  cnew.tmmcInit();
//  mcnew.setNFreqCheckE(10, 1e-7);
//  mcnew.initRestart("tmp/rst", 10);
//  mcnew.nMolSeek(nMolMax, "../forcefield/data.cg3_60_43_1", 1e6);
//  mcnew.runNumTrials(10);
//  p.initEnergy();
//  EXPECT_NEAR(pLJnew.peTot(), p.peTot(), 1e-11);
//
//  p.initEnergy();
//  mc.runNumTrials(10);
//  pLJnew.initEnergy();
//  EXPECT_NEAR(pLJnew.peTot(), p.peTot(), 1e-11);
//
//  // make new wltmmc from restart file, run both 10 steps and check they are the same
//  mcnew.writeRestart("tmp/rst");
//  WLTMMC mcnew2("tmp/rst");
//  EXPECT_EQ(mcnew.space()->nMol(), mcnew2.space()->nMol());
//  EXPECT_NEAR(1., mcnew.pair()->peTot()/mcnew2.pair()->peTot(), 1e-10);
//  mcnew.runNumTrials(10);
//  mcnew2.runNumTrials(10);
//  EXPECT_EQ(mcnew.space()->nMol(), mcnew2.space()->nMol());
//  EXPECT_NEAR(1., mcnew.pair()->peTot()/mcnew2.pair()->peTot(), 1e-10);
//  p.initEnergy();
//}

TEST(MC, b2hardsphere) {
  Space s(3, 0);
  s.addMolInit("../forcefield/data.lj");
  PairHS p(&s, 1);
  CriteriaMetropolis c(1, 1);
  MC mc(&s, &p, &c);
  mc.setNumTrials(1e9);
  mc.initLog("tmp/tmp", 1e5);
  double b2, b2er;
  const double tol = 1e-2;
  mc.b2(tol, b2, b2er);
  EXPECT_NEAR(b2, 2./3.*PI, tol*3);
}

TEST(MC, b2mayer) {
  ranInitByDate();
  Space space(3, 0);
  space.addMolInit("../forcefield/data.lj");
  vector<double> xAdd(space.dimen(), 0.);
  space.xAdd = xAdd;
  space.addMol(0);
  xAdd[0] = 1.1;
  space.xAdd = xAdd;
  space.addMol(0);
  
  PairLJMulti pair(&space, 1e5);
  pair.initData("../forcefield/data.lj");
  pair.cutShift(1);
  pair.initEnergy();
  PairHS pairRef(&space, 1);
  pairRef.initData("../forcefield/data.lj");
  pairRef.initEnergy();

  const double boxl = 2.*(2.*space.maxMolDist() + pair.rCut());
  space.lset(boxl);

  CriteriaMayer crit(0.2);
  crit.initPairRef(&pairRef);
  MC mc(&space, &pair, &crit);
  transformTrial(&mc, "translate", 0.1);
  //mc.setNFreqTune(1e5);
  const int nfreq = 1e0;
  mc.setNFreqCheckE(1, 1e10);
  //mc.setNFreqCheckE(nfreq, 1e-8);
  mc.initLog("tmp/mayer", nfreq);
  mc.initMovie("tmp/mayer", nfreq);
  
  //mc.nMolSeek(2, "../forcefield/data.lj", 1e8);
  
  mc.runNumTrials(1e4);
  //mc.runNumTrials(1e6);
  //cout << "pe " << pair.peTot() << endl;
  EXPECT_NEAR(crit.b2ratio(), 0.2, 1.);
  // beta = 1 takes too long
  // // EXPECT_NEAR(crit.b2ratio(), -2.5381, DTOL);
}

TEST(MC, ljMD) {
  Space space(3,0);
  for (int dim = 0; dim < space.dimen(); ++dim) space.lset(10,dim);
  space.addMolInit("../forcefield/data.atom");
  PairLJ pair(&space, 3.);
  pair.initEnergy();
  CriteriaMetropolis criteria(0.5, 0.1);
  MC mc(&space, &pair, &criteria);

  mc.nMolSeek(50, "../forcefield/data.atom");
  
  shared_ptr<TrialMD> tmd = make_shared<TrialMD>();
  mc.initTrial(tmd);

  mc.runNumTrials(100);
}

TEST(MC, ljWall) {
  feasst::Space space(3,0);
  space.lset(20);
  space.addMolInit("../forcefield/data.lj");
  
  feasst::PairLJMulti pairLJ(&space, 3.);
  pairLJ.initLMPData("../forcefield/data.lj");
  pairLJ.linearShift(1);
  
  feasst::Barrier barrier;
  barrier.addOrthogonalPlanar(5, 1, 0);
  barrier.addOrthogonalPlanar(-5, -1, 0);
  feasst::PairWall pairWall(&space, &barrier);
  
  feasst::PairHybrid pair(&space, space.minl()/2.);
  pair.addPair(&pairLJ);
  pair.addPair(&pairWall);
  pair.initEnergy();
  
  const double beta = 1.;
  feasst::CriteriaMetropolis criteria(beta, exp(-1));
  feasst::MC mc(&space, &pair, &criteria);
  feasst::transformTrial(&mc, "translate", 0.1);
  
  const int nfreq = 1e2;
  mc.initMovie("tmp/wall", nfreq);
  mc.initLog("tmp/wall", nfreq);
  mc.setNFreqCheckE(nfreq, 1e-8);
  mc.setNFreqTune(nfreq);
  mc.nMolSeek(50, "../forcefield/data.lj");

  mc.runNumTrials(1e4);
}

TEST(MC, PairPatchKF) {
  feasst::Space space;
  space.lset(8);
  space.addMolInit("../forcefield/data.onePatch");
  feasst::PairPatchKF pair(&space, 2., 90);
  pair.initEnergy();
  feasst::CriteriaMetropolis criteria(1., exp(-1));
  feasst::MC mc(&space, &pair, &criteria);
  feasst::transformTrial(&mc, "translate", 1);
  feasst::transformTrial(&mc, "rotate", 1);
  space.addTypeForCluster(0);
  feasst::clusterTrial(&mc, "clustertrans", -1, 1.);
  feasst::clusterTrial(&mc, "clusterrotate", -1, 1.);
  feasst::gcaTrial(&mc);
  mc.nMolSeek(12, "../forcefield/data.onePatch");
  pair.initEnergy();
  const int nfreq = 1e0;
  mc.initLog("tmp/patchlog", nfreq); 
  mc.initMovie("tmp/patchmovie", nfreq); 
  mc.setNFreqCheckE(nfreq, feasst::DTOL);
  mc.runNumTrials(1e2);
}

