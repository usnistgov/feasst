import os, sys
feasstdir = os.getenv("FEASST_INSTALL_DIR_") + "/build"
if (not os.path.isfile(feasstdir+"/_feasst.so")):
  feasstdir = os.getenv("FEASST_INSTALL_DIR_") + "/src"
sys.path.append(feasstdir)
import feasst
import math, argparse

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("--openMP", help="use openMP parallelization", dest='openMP', action='store_true')
parser.set_defaults(openMP=False)
parser.add_argument("--npr", "-p", help="number of max trials", default=int(1e10), type=int)
parser.add_argument("--boxl", "-l", help="box length", default=8, type=float)
parser.add_argument("--rCut", "-r", help="cutoff distance", default=3., type=float)
parser.add_argument("--temp", "-t", help="temperature", default=1.2, type=float)
parser.add_argument("--lnz", "-z", help="activity", default=-2.902929, type=float)
parser.add_argument("--nfreq", "-f", help="number of trials per print", default=int(1e5), type=int)
parser.add_argument("--nMolMax", "-x", help="maximum number of mols", default=390, type=int)
parser.add_argument("--molName", "-m", help="molecule file name", default="data.lj", type=str)
args = parser.parse_args()
print(args)

# initialize simulation domain
feasst.ranInitByDate()
space = feasst.Space(3, 0)
for dim in range(space.dimen()): space.lset(args.boxl, dim)
addMolType=space.install_dir()+"/forcefield/"+args.molName
space.addMolInit(addMolType)

# initialize pair-wise interactions
pair = feasst.PairLJ(space, args.rCut)
pair.initEnergy()

# acceptance criteria
nMolMin = 0
criteria = feasst.CriteriaWLTMMC(1./args.temp, math.exp(args.lnz),
  "nmol", nMolMin - 0.5, args.nMolMax + 0.5, args.nMolMax - nMolMin + 1)
criteria.collectInit()
criteria.tmmcInit()

# initialize monte carlo
mc = feasst.WLTMMC(space, pair, criteria)
mc.weight = 3./4.
feasst.transformTrial(mc, "translate")
mc.weight = 1./8.
feasst.deleteTrial(mc)
mc.weight = 1./8.
feasst.addTrial(mc, addMolType)

# if using parallelization, allow configuration swaps between processors
if args.openMP:
  mc.weight = 1./args.nMolMax
  mc.confSwapTrial()

# output log, lnpi and movie
mc.initLog("log", args.nfreq)
mc.initColMat("colMat", args.nfreq)
mc.setNFreqCheckE(args.nfreq, 1e-8)
mc.setNFreqTune(args.nfreq)
mc.initMovie("movie", args.nfreq)
#mc.initXTC("movie", args.nfreq)
mc.initRestart("tmp/rst", args.nfreq)

# production tmmc simulation
if not args.openMP:
  mc.runNumTrials(args.npr)
else:
  mc.initWindows(2.,  # exponent that determines size of windows
                 0)   # extra macrostate overlap between processors
  mc.runNumSweeps(1,  # number of "sweeps"
                 -1)  # maximum number of trials. Infinite if "-1".

