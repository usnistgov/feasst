import pandas as pd
import argparse
parser = argparse.ArgumentParser()
required = parser.add_argument_group('required arguments')
required.add_argument("--file_name", "-f", help="file name", type=str, required=True)
args = parser.parse_args()

print(args.file_name)

# find the line which begins with the characters 'macrostate'
with open(args.file_name, 'r') as f:
    lines = f.readlines()
    count = 0
    for line in lines:
        if line.split(',')[0] == 'macrostate':
            header = count
        count += 1

df = pd.read_csv(args.file_name, header=header)
#print(df)

import matplotlib.pyplot as plt
plt.plot(df['macrostate'], df['ln_prob'])
plt.show()
