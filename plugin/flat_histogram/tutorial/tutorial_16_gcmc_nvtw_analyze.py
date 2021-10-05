import pandas as pd
import feasst as fst

min_particles=0
max_particles=370
data=list()
for num in range(min_particles, max_particles+1):
#for num in reversed(range(min_particles, max_particles+1)):
    df = pd.read_csv('crit'+str(num)+'.txt', header=1)
    print(df['c0'].values[0])
    #print(df['visits'])
    #print(df['c0'], df['c1'], df['c2'])
    data.append([[float(df['c0'].values[0]),
                  float(df['c1'].values[0]),
                  float(df['c2'].values[0])]])

print(data)
tm = fst.TripleBandedCollectionMatrix(data)
lnp = fst.LnProbability()
lnp.resize(len(data))
tm.compute_ln_prob(lnp)
print(lnp.values())

import matplotlib.pyplot as plt
plt.plot(lnp.values())
plt.show()

