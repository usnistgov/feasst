import pandas as pd
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
required = parser.add_argument_group('required arguments')
required.add_argument("--file_name", "-f", help="file name", action='append', required=True)
args = parser.parse_args()
print(args.file_name)

# find the line which begins with the characters 'state'
for fname in args.file_name:
    with open(fname, 'r') as f:
        lines = f.readlines()
        count = 0
        for line in lines:
            if line.split(',')[0] == 'state':
                header = count
            count += 1
    df = pd.read_csv(fname, header=header)
    #print(df)
    plt.plot(df['state'], df['ln_prob'])
plt.legend()
plt.show()
