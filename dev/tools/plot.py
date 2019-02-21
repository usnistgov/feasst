''' use this script to visualize test of random unit sphere '''

import pandas as pd
df = pd.read_csv("tt", delim_whitespace=True, header=None)
print(df)
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(df[0], df[1], df[2])

plt.show()
