"""@docstring
@brief Tests pair interactions
@author Nathan A. Mahynski and Harold Hatch
@date 5/16/2017
@filename test_pair.py
"""

import sys, os
import unittest
import numpy as np

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
