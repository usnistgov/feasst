import pandas as pd
import matplotlib.pyplot as plt

#df = pd.read_csv("wl/crit.txt", delim_whitespace=True, header=3)
#plt.plot(df["macrostate"], df["ln_prob"], label="wl")
df = pd.read_csv("crit.txt", delim_whitespace=True, header=3)
plt.plot(df["macrostate"], df["ln_prob"])
#plt.plot(df["macrostate"], df["c0"])
#plt.plot(df["macrostate"], df["c1"])
#plt.plot(df["macrostate"], df["c2"])
plt.legend()
plt.show()
