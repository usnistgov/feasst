#! /usr/bin/env python

import os, sys
#from libxdrfile import xdrfile_open, xdrfile_close, read_xtc_natoms, read_xtc, DIM, exdrOK
feasstdir = os.getenv("HOME") + "/jeetain/feasstv0.1.3/build"
sys.path.append(feasstdir)
import feasst
import math, argparse, time
start_time = time.time()

# parse arguments and print to log file
parser = argparse.ArgumentParser()
parser.add_argument("--npr", "-r", help="number of max trials", default=int(1e6), type=int)
parser.add_argument("--nPart", "-n", help="number of particles", default=38, type=int)
parser.add_argument("--pres", "-p", help="pressure", default=1., type=float)
parser.add_argument("--za", "-a", help="activity of mol type A", default=math.exp(-4.568214), type=float)
parser.add_argument("--zb", "-b", help="activity of mol type B", default=math.exp(-3), type=float)
parser.add_argument("--boxl", "-l", help="box length", default=8, type=float)
parser.add_argument("--temp", "-t", help="temperature", default=1.5, type=float)
parser.add_argument("--nfreq", "-f", help="number of trials per output", default=int(1e4), type=int)
parser.add_argument("--molNameA", "-x", help="molecule file name", default="data.ljA", type=str)
parser.add_argument("--molNameB", "-y", help="molecule file name", default="data.ljB", type=str)
args = parser.parse_args()
print(args)

# initialize random number seed
feasst.ranInitByDate()

# initialize simulation domain
space = feasst.Space(3, 0)
for dim in range(space.dimen()): space.lset(args.boxl, dim)

# initialize pair-wise interactions
rCut = 3
space.updateCells(rCut)
pair = feasst.PairLJMulti(space, rCut)
space.addMolInit(args.molNameA)
pair.initLMPData(args.molNameA)
space.addMolInit(args.molNameB)
pair.initLMPData(args.molNameB)
pair.setLambdaij(0,0,1) #lambda_AA
pair.setLambdaij(0,1,0.5) #lambda_AB
pair.setLambdaij(1,1,1) #lambda_BB
#pair.epsijset(0, 0, 0.95)
#pair.epsijset(1, 1, 0.95)
pair.linearShift(1)
#pair.cutShift(1)
pair.initEnergy()

# acceptance criteria
criteria = feasst.CriteriaMetropolis(1./args.temp, args.za)
criteria.addActivity(args.zb)
criteria.pressureset(args.pres)

# initialize monte carlo
mc = feasst.MC(space, pair, criteria)

# translation trial move
mc.weight = 3./4.;
maxMoveParam = 0.1
feasst.transformTrial(mc, "translate", maxMoveParam)

# generate initial configuration
assert(args.nPart % 2 == 0)
mc.nMolSeek(args.nPart/2, args.molNameA)
mc.nMolSeek(args.nPart, args.molNameB)

## volume change trial
mc.weight = 1./args.nPart;
maxMoveParam = 0.001
feasst.transformTrial(mc, "vol", maxMoveParam)

# identity swap trial
mc.weight = 1./4.;
trialSwap = feasst.TrialSwap(args.molNameA, args.molNameB)
mc.initTrial(trialSwap)

# output log and movie
mc.initLog("log", args.nfreq)
mc.setNFreqCheckE(args.nfreq, 1e-8)
mc.setNFreqTune(args.nfreq)
mc.initMovie("movie", args.nfreq) 
mc.initXTC("moviex", args.nfreq) 
mc.initRestart("tmp/rst", args.nfreq)
  
# production simulation
mc.initProduction()
xA = feasst.Accumulator()
trialPerAccum = 1
for step in range(args.npr/trialPerAccum):
  mc.runNumTrials(trialPerAccum)
  xA.accumulate(float(space.nMolType()[0])/
                float(space.nMol()))

print("xA", xA.average(), "+/-", xA.blockStdev())

