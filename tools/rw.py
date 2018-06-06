# Reweight a collection matrix file to a given lnz, if given.
# If a lnz is not given, then attempt to find saturation.
#   For saturation, if a volume is not given, attempt to find the volume from a
#   restart file in order to compute the density and pressure at saturation.

import os
import math
import feasst
import argparse

LNZDEFAULT=11234533

parser = argparse.ArgumentParser()
parser.add_argument('--inFile',  '-i', help="input collection matrix file", default="colMat.txt",   type=str)
parser.add_argument('--outFile', '-o', help="output collection matrix file",default="colMatrw.txt", type=str)
parser.add_argument('--lnz',     '-z', help="ln(activity)",                 default=LNZDEFAULT,       type=float)
parser.add_argument('--phaseBoundary', '-p', help="assign number of molecules as phase boundary", default=-1, type=int)
parser.add_argument('--volume', '-v', help="volume of simulation domain", default=-1, type=float)
parser.add_argument('--rstFile', '-r', help="restart file for space domain", default="tmp/rstspace", type=str)
args = parser.parse_args()
print(args)

# if the lnz was set to some value, reweight to the given lnz
if (args.lnz != LNZDEFAULT):
  criteria = feasst.CriteriaWLTMMC(args.inFile)
  criteria.readCollectMat(args.inFile)
  criteria.lnPIrw(math.exp(args.lnz))
  criteria.printRWinit();
  criteria.printCollectMat(args.outFile)

# otherwise, attempt to find saturation and reweight the lnzsat
else:
  # if volume is not given or unphysical, attempt to find volume from restart
  if (args.volume < 0):
    assert(os.path.isfile(args.rstFile))  # provide either -v or -r flags
    space = feasst.Space(args.rstFile)
  else:
    space = feasst.Space()
    space.initBoxLength(1)
    space.scaleDomain(args.volume)

  pair = feasst.PairIdeal(space)
  criteria = feasst.CriteriaWLTMMC(args.inFile)
  wltmmc = feasst.WLTMMC(space, pair, criteria)
  criteria.readCollectMat(args.inFile)
  if (args.phaseBoundary != -1): criteria.setPhaseBoundary(args.phaseBoundary)
  wltmmc.printSat();


