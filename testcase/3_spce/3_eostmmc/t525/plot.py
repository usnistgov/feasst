import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("colMat", delim_whitespace=True, comment="#", header=None)
print(df.head())

plt.plot(df[0], df[1])
plt.show()

