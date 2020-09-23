"""Obtain average and standard deviation from the lj.msdb.* files"""

import os
import pandas as pd
import numpy as np

temps = ["300"]
for temp in temps:
    os.system('tar -xf c2.' + str(temp) + 'K.tar.gz')
    dfl=list()
    for replica in range(4):
        df=pd.read_csv("c2."+temp+"K.r"+str(replica+1)+".p_macro.dat", delim_whitespace=True, header=None)
        # df["replica"] = replica
        df.rename(columns={0: "N", 1: "lnPI"}, inplace=True)

        dfl.append(df)

    df=pd.concat(dfl)
    #df.to_csv("agg.csv")
    grp=df.groupby("N")
    dfstd = grp.std()
    dfstd.rename(columns={"energy": "energystd", "lnPI": "lnPIstd"}, inplace=True)
    df=pd.concat([grp.mean(), dfstd], axis=1)
    df['lnPI'] -= np.log(sum(np.exp(df['lnPI'])))
    df.to_csv("stat_c2_" + str(temp) + ".csv")
