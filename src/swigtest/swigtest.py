#! /usr/bin/env python
import sys, os
feasstdir = "../build"
if (not os.path.isfile(feasstdir+"/_feasst.so")):
  feasstdir = ""
print(feasstdir)
sys.path.append(feasstdir)
import feasst
import unittest
import numpy as np

class TestFunctions(unittest.TestCase):
  #check that random number generator selects particles randomly
  def testMyRanNum(self):
    imin=5
    imax=9
    n=10000
    x=np.array(range(imax+1))
    feasst.ranInitByDate()
    s=feasst.Space(3,0)
    for i in range(n):
        ran=s.uniformRanNum(imin, imax)
        x[ran]=x[ran]+1
    for i in range(imin,imax+1):
        self.assertAlmostEqual(1./(imax-imin+1),x[i]/float(n), 1)

class TestSpace(unittest.TestCase):
    #check that space class returns the correct number of atoms
    def testNatoms(self):
        for natoms in range(1,10):
            for dimen in range(1,4):
                s=feasst.Space(dimen,0)
                self.assertEqual(dimen,s.dimen())

    #check that init_config works as expected
    def testInitConfandXaccessor(self):
        for natoms in range(1,10):
            for dimen in range(1,4):
                s=feasst.Space(dimen,0)
                s.init_config(natoms)
                self.assertEqual((natoms-1)*dimen*0.95,s.x(natoms-1,0))

    #checks setting positions of particles
    def testXset(self):
        s=feasst.Space(3,0)
        s.init_config(3)
        x0=4.
        s.xset(x0,0,0)
        self.assertEqual(x0,s.x(0,0))
        
    #checks setting box dimensions
    def testLset(self):
        s=feasst.Space(3,0)
        l0=4.
        self.assertEqual(0.,s.l(0))
        s.lset(4.,0)
        self.assertEqual(l0,s.l(0))
        

class TestPair(unittest.TestCase):
    #check that pair class returns the correct potential cut-off
    def testRcut(self):
        rCut=1
        s=feasst.Space(2,0)
        s.init_config(2)
        p=feasst.PairLJ(s,rCut)
        self.assertEqual(rCut,p.rCut())

#    #check LJ energy and virial of configurations from standard reference sim web
#    # http://www.nist.gov/mml/csd/informatics_research/lj_refcalcs.cfm
#    def testSRSWLJ(self):
#      for conf in range(1,5):
#        for lrc in [0, 1]:
#          rCut=3
#          if (conf==1):
#              boxl=10
#              peSRSW=-4351.5
#              peSRSWfull=-4549.99
#              peplaces=1
#              vrSRSW=-568.67
#              vrplaces=2
#          if (conf==2):
#              boxl=8
#              peSRSW=-690.0
#              peSRSWfull=-714.23
#              peplaces=1
#              vrSRSW=-568.46
#              vrplaces=2
#          if (conf==3):
#              boxl=10
#              peSRSW=-1146.7
#              peSRSWfull=-1196.322
#              peplaces=1
#              vrSRSW=-1164.9
#              vrplaces=1
#          if (conf==4):
#              boxl=8
#              peSRSW=-16.790
#              peSRSWfull=-17.33517
#              peplaces=3
#              vrSRSW=-46.249
#              vrplaces=3
#          s=feasst.Space(3,0)
#          s.addMolInit("forcefield/data.lj")
#          s.lset(boxl)
#          s.readxyz("test/lj/srsw/lj_sample_config_periodic"+str(conf)+".xyz")
#          p=feasst.PairLJ(s,rCut)
#          p.initLMPData("forcefield/data.lj")
#          p.lrcFlag=lrc
#          p.Forces()
#          print "conf", conf, "pe", p.peTot(), "N", s.nMol(), "test/lj/srsw/lj_sample_config_periodic"+str(conf)+".xyz", "natom", s.natom()
#          if (lrc==0):
#            self.assertAlmostEqual(peSRSW, p.peTot(), peplaces)
#            self.assertAlmostEqual(vrSRSW, p.vrTot(), vrplaces)
#          if (lrc==1):
#            self.assertAlmostEqual(peSRSWfull, p.peTot(), peplaces)
#
#    #check SPCE energy of configurations from standard reference sim web
#    # http://www.nist.gov/mml/csd/informatics_research/lj_refcalcs.cfm
#    def testSRSWSPCE(self):
#      for conf in range(1,2):
#        for lrc in [0, 1]:
#          rCut=10
#          if (conf==1):
#              peplaces=6
#              boxl=20
#              peSRSWLJ=9.95382E+04
#              peSRSWLRC=-8.23710E+02
#              peSRSWReal=-5.58889E+06
#              peSRSWFrr=6.27797E+03
#              peSRSWSelf=-2.84469E+06+2.80999E+06
#              peSRSW=-4.88596E+05
#          s=feasst.Space(3,0)
#          for i in range(s.dimen()):
#              s.lset(boxl, i)
#          s.readxyz("test/spce/srsw/spce_sample_config_periodic"+str(conf)+".xyz")
#          p=feasst.PairLJ(s,rCut)
#          p.lrcFlag=lrc
#          p.Forces()
#          if (lrc==0):
#            self.assertAlmostEqual(peSRSW, p.peTot(), peplaces)
#            self.assertAlmostEqual(vrSRSW, p.vrTot(), vrplaces)
#          if (lrc==1):
#            self.assertAlmostEqual(peSRSWfull, p.peTot(), peplaces)

class TestTrial(unittest.TestCase):
#    #check wltmmc algorithm with ideal gas
#    def testTrialIdealTMMC(self):
#        feasst.ranInitByDate()
#        s=feasst.Space(3,0)
#        s.init_config(12)
#        p=feasst.PairIdeal(s,5)
#        rAbove=3
#        rBelow=1
#        p.initNeighList(rAbove, rBelow)
#        p.buildNeighList()
#        #c=feasst.CriteriaMetropolis(1,0.0000001)
#        c=feasst.CriteriaWLTMMC(1,0.0000001,"nmol",0-0.5,s.nMol()+0.5,s.nMol()+1)
#        tt=feasst.TrialTransform(s,p,c,"translate")
#        ta=feasst.TrialAdd(s,p,c,"lj")
#        td=feasst.TrialDelete(s,p,c)
#        tavb1=feasst.TrialAVB(s,p,c, 0.9, rAbove, rBelow, 1)
#        tavb2=feasst.TrialAVB(s,p,c, 0.9, rAbove, rBelow, 2)
#        tavb3=feasst.TrialAVB(s,p,c, 0.9, rAbove, rBelow, 3)
#        tt.maxMoveParam = 5
#        p.Forces()
#        for i in range(0,30000):
#            tt.attempt()
#            ta.attempt()
#            td.attempt()
#            tavb1.attempt()
#            tavb2.attempt()
#            tavb3.attempt()
#            #print i, s.natom()
#            print i, s.natom(), c.lnf()

    #check that for all pairs and criteria, neighbor list and potential energy are updated properly
    def testTrialPEandNeighUpdate(self):
        feasst.ranInitByDate()
        ctlist = ["metropolis","wltmmc"]
        ptlist = ["lj","ideal"]
        for pt in ptlist:
          for ct in ctlist:
            print(pt, ct)
            s=feasst.Space(3,0)
            s.init_config(12)
            s.addMolInit("../forcefield/data.atom")
            if (pt == "ideal"):
              p=feasst.PairIdeal(s,5)
            elif (pt == "lj"):
              p=feasst.PairLJ(s,5)
            rAbove=3
            rBelow=1
            p.initNeighList(rAbove, rBelow)
            p.buildNeighList()
            if (ct == "metropolis"):
              c=feasst.CriteriaMetropolis(1,0.01)
            elif (ct == "wltmmc"):
              c=feasst.CriteriaWLTMMC(1,0.01,"nmol",0-0.5,s.nMol()+0.5,s.nMol()+1)
            mc=feasst.MC(s, p, c)
            feasst.transformTrial(mc, "translate", 5)
            feasst.deleteTrial(mc)
            feasst.addTrial(mc, "../forcefield/data.atom")
            #mc.avbTrial(0.9, rAbove, rBelow, 1)
            #mc.avbTrial(0.9, rAbove, rBelow, 2)
            #mc.avbTrial(0.9, rAbove, rBelow, 3)
            p.initEnergy()
            mc.runNumTrials(300)
            peTot=p.peTot()
            #peLJ=p.peLJ()
            p.Forces()   
            self.assertAlmostEqual(p.peTot(), peTot, 12)
            #self.assertAlmostEqual(p.peLJ(), peLJ, 13)
            self.assertEqual(1, p.checkNeigh())

class TestMC(unittest.TestCase):
  def testMC(self):
    s=feasst.Space(3,0)
    s.lset(10)
    p=feasst.PairLJMulti(s, 3)
    c=feasst.CriteriaWLTMMC(1, 0.1, "nmol",0.5,1.5,2)
    mc=feasst.WLTMMC(s, p, c)

if __name__ == '__main__':
    unittest.main()


