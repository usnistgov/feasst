import pandas as pd
import matplotlib.pyplot as plt

np=1

dfall = list()
prev=0
for p in range(np):
  df1=pd.read_csv("crit"+str(p)+".txt", delim_whitespace=True, header=3)
  df1["ln_prob"] += prev - df1["ln_prob"].iloc[0]
  dfall.append(df1)
  prev = df1["ln_prob"].iloc[-1]

df = pd.concat(dfall)
df.to_csv("combined.csv")

plt.plot(df["macrostate"], df["ln_prob"])
plt.legend()
plt.show()
