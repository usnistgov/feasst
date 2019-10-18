"""Create a plot from a criteria file"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--prepend", "-f", help="criteria prepend name", default="crit", type=str)
parser.add_argument("--num_processors", help="number of processors", default=1, type=int)
parser.add_argument("--label", help="analyze label", default="ln_prob", type=str)
args = parser.parse_args()
print("#", args)

def combine_criteria():
    """Combine multiple criteria files together, assuming that the last macrostate in the ith
    file is equivalent to the first macrostate in the (i+1)th file.
    """
    dfall = list()
    prev = 0
    for proc in range(args.num_processors):
        df1 = pd.read_csv(args.prepend+str(proc)+".txt", delim_whitespace=True, header=3)
        df1[args.label + "_combine"] = df1[args.label]
        df1[args.label + "_combine"] += prev - df1[args.label].iloc[0]
        dfall.append(df1)
        prev = df1[args.label].iloc[-1]
    return pd.concat(dfall)

df = combine_criteria()
df.to_csv("combined.csv")
plt.scatter(df["macrostate"], df[args.label + "_combine"])
plt.xlabel("macrostate", fontsize=16)
ylabel = args.label
if args.label == 'ln_prob':
    ylabel = r'$\ln\Pi$'
plt.ylabel(ylabel, fontsize=16)
if args.num_processors > 1:
    plt.legend()
plt.savefig("crit.svg")
plt.show()
