"""Create a plot of energy from state files"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import feasst
import pyfeasst

parser = argparse.ArgumentParser()
parser.add_argument("--prepend", "-f", help="state file prepend name", default="state", type=str)
parser.add_argument("--num_processors", help="number of processors", default=1, type=int)
parser.add_argument("--label", help="analyze label", default="energy", type=str)
args = parser.parse_args()
print("#", args)

def combine_energy():
    """Aggregate individual files into one data frame"""
    dfall = list()
    for proc in range(args.num_processors):
        monte_carlo = feasst.MonteCarlo().deserialize(
            pyfeasst.read_checkpoint("checkpoint" + str(proc) + ".txt"))
        for state in range(monte_carlo.criteria().num_states()):
            df = pd.read_csv(args.prepend + str(state) + "_energy" + str(proc) + ".txt", delim_whitespace=True)
            df['state'] = state
            dfall.append(df)
    return pd.concat(dfall)

en = combine_energy()
dfsrsw = pd.read_csv("../test/data/stat150.csv")

plt.xlim([-1, 6])
plt.ylim([-0.4, 0.1])
plt.errorbar(en["state"], en["average"], yerr=en["block_stdev"], fmt='o', label="FEASST")
plt.errorbar(dfsrsw["N"], dfsrsw[args.label], yerr=dfsrsw[args.label + "std"], fmt='x', label="NIST SRSW")
plt.xlabel("macrostate", fontsize=16)
plt.ylabel(args.label, fontsize=16)
plt.legend()
plt.savefig(args.label + ".svg")
plt.show()
