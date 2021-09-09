import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("~/feasst/build/tmp/butane_dihedral.txt")

def en(phi):
    c0=0
    c1=2.951883663322944 #355.03
    c2=-0.5669632059318694 #-68.19
    c3=6.579400558997022 #791.32
    return c0+c1*(1+np.cos(phi))+c2*(1-np.cos(2*phi))+c3*(1+np.cos(3*phi))

plt.plot(df['bin'], en(df['bin']))
plt.show()

print(df['bin'].values)
print(np.cos(df['bin'].values))
prob = np.exp(-en(df['bin'].values))
prob /= np.sum(prob)
df['num'] /= np.sum(df['num'])

plt.plot(df['bin'], df['num'])
plt.plot(df['bin'], prob)

plt.show()
