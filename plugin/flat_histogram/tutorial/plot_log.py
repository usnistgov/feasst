"""Create a plot from the log file."""

import argparse
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--filename", "-f", help="log file name", default="log0.txt", type=str)
parser.add_argument("-x", help="column name for abscissa", default="attempt", type=str)
parser.add_argument("-y", help="column name for ordinate", default="state", type=str)
args = parser.parse_args()
print("#", args)

df = pd.read_csv(args.filename, delim_whitespace=True, header=0)
plt.scatter(df[args.x], df[args.y])
plt.xlabel(args.x, fontsize=16)
plt.ylabel(args.y, fontsize=16)
plt.savefig(args.filename + ".svg")
plt.show()
