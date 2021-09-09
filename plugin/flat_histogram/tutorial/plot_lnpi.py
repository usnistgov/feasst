import pandas as pd
import argparse
import matplotlib.pyplot as plt
import feasst as fst
import pyfeasst

parser = argparse.ArgumentParser()
parser.add_argument("--file_name", "-f", type=str, help="file name", action='append', required=True)
parser.add_argument("--label", "-l", type=str, help="variable label", default='ln_prob')
parser.add_argument("--delta_conjugate", "-d", type=float, help="reweight by change in conjugate variable", default=0.)
args = parser.parse_args()
print(args.file_name)

for fname in args.file_name:
    df = pd.read_csv(fname, header=pyfeasst.line_beginning_with_state(fname))
    gce = fst.GrandCanonicalEnsemble(
        fst.Histogram(fst.args({"width": "1", "max": str(df["state"].iloc[-1]),
                                              "min": str(df["state"].iloc[0])})),
        fst.LnProbability(fst.DoubleVector(df["ln_prob"])))
    #print(df)
    #plt.plot(df['state'], df[args.label], label=fname)
    gce.reweight(args.delta_conjugate)
    plt.plot(df['state'], gce.ln_prob().values(), label=fname)
plt.legend()
plt.show()
