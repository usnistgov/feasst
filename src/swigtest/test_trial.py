"""
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of this agreement (see LICENSE.txt) and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
"""

"""@docstring
@brief Tests trial moves
@author Nathan A. Mahynski
@date 5/16/2017
@filename test_trail.py
"""

import sys, os
import unittest
import numpy as np

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
				s.addMolInit(os.path.dirname(os.path.realpath(__file__))+"/../../forcefield/data.atom")
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
				feasst.addTrial(mc, os.path.dirname(os.path.realpath(__file__))+"/../../forcefield/data.atom")
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

if __name__ == "__main__":
	def find_all(name, path):
	    result = []
	    for root, dirs, files in os.walk(path):
	        if name in files:
	            result.append(os.path.join(root, name))
	    return result

	feasstname = find_all("_feasst.so", os.path.dirname(os.path.realpath(__file__))+"/../../")
	if (len(feasstname) > 1):
		print ('Found multiple FEASST sources which may cause an error: ', feasstname)
		sys.exit()
	else:
		feasstdir = feasstname[0].split('/_feasst.so')[0]
		sys.path.append(feasstdir)

	import feasst
	unittest.main()
