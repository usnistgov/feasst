import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--num_procs", "-n", help="number of processors", type=int, required=True)
parser.add_argument("--prepend", help="file name begins with these characters (e.g., \"en0.txt\" has a prepend of \"en\"", type=str, default="en")
parser.add_argument("--append", help="file name ends with these characters (e.g., \"en0.txt\" has an append of \".txt\"", type=str, default=".txt")
args = parser.parse_args()

dr=''
ens=list()
for proc in range(args.num_procs):
    file1 = open(dr+'crit'+str(proc)+'.txt', 'r')
    Lines = file1.readlines()
    exec('params={'+Lines[0][1:]+'}')
    #print(params)
    #print(params['soft_min'])
    en = pd.read_csv('en'+str(proc)+'.txt')
    mn=params['soft_min']
    if proc == 0:
        mn = 0
    mx=params['soft_max']
    if proc == args.num_procs - 1:
        mx = en['state'].values[-1]
#    print(en[mn:mx+1])
    ens.append(en[mn:mx+1])

df=pd.concat(ens)
#print(df)
df.to_csv('en.txt')

#plt.errorbar(df['state'], df['average'], df['block_stdev'])
#
#srsw=pd.read_csv('~/feasst/plugin/flat_histogram/test/data/stat150.csv')
#print(srsw)
#plt.errorbar(srsw['N'], srsw['energy'], srsw['energystd'])
#
#plt.show()

#nconfig=int(len(Lines)/373)
#for config in range(nconfig):
#    exec('params={'+Lines[config*373][1:]+'}')
#df=pd.read_csv(args.prepend+"0"+args.append)
#for proc in range(1, args.num_procs):
#    df2=pd.read_csv(args.prepend+str(proc)+args.append)
#    df = df+df2
#    #df.append(df)
#print(df['moment0'])
#print(df['moment1'])
#df['average'] = df['moment1']/df['moment0']
#print(df['average'])
##print(df)
#
