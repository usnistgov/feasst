# Obtain average and standard deviation from the lj.msdb.* files
import pandas as pd
dfl=list()
temp="150"
for replica in range(3):
  # read energy
  df=pd.read_csv("lj.msdb.t"+temp+"."+str(replica+1)+".energy.dat", delim_whitespace=True, header=None)
  df["replica"] = replica
  df.rename(columns={0: "N", 1: "energy"}, inplace=True)

  # add macro
  dfmacro=pd.read_csv("lj.msdb.t"+temp+"."+str(replica+1)+".p_macro.dat", delim_whitespace=True, header=None)
  df["lnPI"] = dfmacro[1]

  dfl.append(df)

df=pd.concat(dfl)
df.to_csv("agg.csv")
grp=df.groupby("N")
df=pd.concat([grp.mean(), grp.std()], axis=1)
df.to_csv("stat.csv")

