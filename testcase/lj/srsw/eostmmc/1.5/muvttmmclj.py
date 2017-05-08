#! /usr/bin/env python

import os, sys
#from libxdrfile import xdrfile_open, xdrfile_close, read_xtc_natoms, read_xtc, DIM, exdrOK
feasstdir = os.getenv("FEASST_INSTALL_DIR_") + "/build"
if (not os.path.isfile(feasstdir+"/_feasst.so")):
  feasstdir = os.getenv("FEASST_INSTALL_DIR_") + "/src"
sys.path.append(feasstdir)
import feasst
import math, argparse

# parse arguments and print to log file
parser = argparse.ArgumentParser()
parser.add_argument("--npr", "-p", help="number of max trials", default=int(1e10), type=int)
parser.add_argument("--boxl", "-l", help="box length", default=8, type=float)
parser.add_argument("--rCut", "-r", help="cutoff distance", default=3., type=float)
parser.add_argument("--temp", "-t", help="temperature", default=1.5, type=float)
parser.add_argument("--lnz", "-z", help="activity", default=-1.568214, type=float)
parser.add_argument("--nfreq", "-f", help="number of trials per print", default=int(1e5), type=int)
parser.add_argument("--nMolMax", "-x", help="maximum number of mols", default=370, type=int)
parser.add_argument("--molName", "-m", help="molecule file name", default="data.lj", type=str)
args = parser.parse_args()
name = os.path.splitext(os.path.basename(__file__))[0]
file = open(name+".log",'w')
print >>file, args

# initialize simulation domain
feasst.ranInitByDate()
space = feasst.Space(3, 0)
for dim in range(space.dimen()): space.lset(args.boxl, dim)
space.addMolInit(args.molName)

# initialize pair-wise interactions
pair = feasst.PairLJ(space, args.rCut)
pair.initEnergy()

# acceptance criteria
nMolMin = 0
criteria = feasst.CriteriaWLTMMC(1./args.temp, math.exp(args.lnz),"nmol",nMolMin-0.5,args.nMolMax+0.5,args.nMolMax-nMolMin+1)

# initialize monte carlo
mc = feasst.WLTMMC(space, pair, criteria)
mc.weight = 3./4.
#tt = feasst.TrialTransform(space, pair, criteria, "translate")
#tt = feasst.TrialTransform("translate")
#mc.initTrial(tt)
#mc.initTrial(feasst.TrialTransform("translate"))
feasst.transformTrial(mc, "translate")
#mc.transformTrial("translate")
mc.weight = 1./8.
#td = feasst.TrialDelete()
#mc.initTrial(td)
#mc.initTrial(feasst.TrialDelete())
feasst.deleteTrial(mc);
#ta = feasst.TrialAdd(args.molName)
#mc.initTrial(ta)
#mc.initTrial(feasst.TrialAdd(args.molName))
feasst.addTrial(mc, args.molName)

# output log, lnpi and movie
mc.initLog("log", args.nfreq)
mc.initColMat("colMat", args.nfreq) 
mc.setNFreqCheckE(args.nfreq, 1e-8)
mc.setNFreqTune(args.nfreq)
mc.initMovie("movie", args.nfreq) 
#mc.initXTC("movie", args.nfreq) 
mc.initRestart("tmp/rst", args.nfreq)
  
# production tmmc simulation
criteria.collectInit()
criteria.tmmcInit()
mc.runNumTrials(args.npr)
#mc.initWindows(1)
#mc.runNumSweeps(1, -1)

