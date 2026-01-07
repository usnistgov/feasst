import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import pandas as pd
df = pd.read_csv('tmp/gauss.txt', header=None)
print(df)
fig = plt.figure()
axes = fig.add_subplot(111, projection='3d')
axes.plot_trisurf(df[0], df[1], df[3])
plt.show()
