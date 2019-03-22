import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("log.txt", delim_whitespace=True, header=None, comment="#")
plt.plot(df[1], df[5])
#plt.plot(df["attempt"], df["N0"])
plt.show()
