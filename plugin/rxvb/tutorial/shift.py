L=24.555
hl=L/2.
import pandas as pd
df = pd.read_csv('tt', sep='\s+', header=None)
df[2] -= hl
df[3] -= hl
df[4] -= hl
df.to_csv('tt2', sep=' ', index=False)
