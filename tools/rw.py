#!/usr/bin/env python

import math, os, sys
feasstdir = os.getenv("FEASST_INSTALL_DIR_") + "/build"
if (not os.path.isfile(feasstdir+"/_feasst.so")):
  feasstdir = os.getenv("FEASST_INSTALL_DIR_") + "/src"
sys.path.append(feasstdir)
import feasst
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--inFile',  '-i', help="input collection matrix file", default="colMat.txt",   type=str)
parser.add_argument('--outFile', '-o', help="output collection matrix file",default="colMatrw.txt", type=str)
parser.add_argument('--lnz',     '-z', help="ln(activity)",                 default=11234533,       type=float)
parser.add_argument('--phaseBoundary', '-p', help="assign number of molecules as phase boundary", default=-1, type=int)
args = parser.parse_args()
print args

if (os.path.isfile("tmp/rstspace")):
  space = feasst.Space("tmp/rstspace")
else:
  space = feasst.Space()
pair = feasst.PairIdeal(space, 0)
criteria = feasst.CriteriaWLTMMC(args.inFile)
wltmmc = feasst.WLTMMC(space, pair, criteria)
criteria.readCollectMat(args.inFile)

if (args.phaseBoundary != -1): criteria.setPhaseBoundary(args.phaseBoundary)
if (args.lnz != 11234533):
  criteria.lnPIrw(math.exp(args.lnz))
  criteria.printRWinit();
  criteria.printCollectMat(args.outFile)
else:
  wltmmc.printSat();


