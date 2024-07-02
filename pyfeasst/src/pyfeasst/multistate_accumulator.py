"""
This module analyzes and manipulates multistate accumulators.
For example, the canonical ensemble average energy for a given number of particles.
"""

from pathlib import Path
import pandas as pd

def splice_by_max_column(prefix, suffix, column='moment0', crit_prefix=None, crit_suffix=None):
    """
    Combine all files with matching prefix and suffix, each with the same number of states,
    by assigning rows based on the maximum in the given column.

    :param str prefix: Find all files beginning with this prefix.
    :param str suffix: Find all files endding with this suffix.
    :param str column:
        Use the row with the maximum column. For example, moment0 is the highest number of samples.
    :param str crit_prefix: Use criteria to read only specific macrostates, if given.
    :param str crit_suffix: Use criteria to read only specific macrostates, if given.

    >>> import pandas as pd
    >>> from pyfeasst import multistate_accumulator
    >>> spliced = multistate_accumulator.splice_by_max_column(prefix="../../tests/lj_enn0s",
    ...                                                       suffix='.txt')
    >>> en0 = pd.read_csv('../../tests/lj_enn0s00.txt')
    >>> en1 = pd.read_csv('../../tests/lj_enn0s01.txt')
    >>> assert spliced['average'][136] == en0['average'][136]
    >>> assert spliced['average'][136] != en1['average'][136]
    >>> assert spliced['average'][137] == en0['average'][137]
    >>> assert spliced['average'][137] == en1['average'][137]
    >>> assert spliced['average'][138] != en0['average'][138]
    >>> assert spliced['average'][138] == en1['average'][138]
    """
    spliced = None
    if crit_prefix is None:
        for filename in Path('.').rglob(prefix+'*'+suffix):
            if spliced is None:
                spliced = pd.read_csv(filename)
            else:
                addition = pd.read_csv(filename)
                newer = addition[addition[column] > spliced[column]]
                if len(newer) > 0:
                    spliced.loc[newer['state'].values[0]: newer['state'].values[-1]] = newer
    else:
        rows = list()
        for filename in sorted(Path('.').rglob(crit_prefix+'*'+crit_suffix)):
            print('filename', filename)
            with open(filename, 'r') as file1:
                lines = file1.readlines()
            exec('iprm={' + lines[0][1:] + '}', globals())
            rows.append([iprm['soft_min'], iprm['soft_max']])
        #print('rows', rows)
        ifile = 0
        lines=list()
        for filename in sorted(Path('.').rglob(prefix+'*'+suffix)):
            #print('filename', filename, 'ifile', ifile)
            with open(filename, 'r') as file1:
                tlines = file1.readlines()
                if ifile == 0: lines += tlines[0]
                lines += tlines[rows[ifile][0]+1:rows[ifile][1]+2]
            ifile += 1
        #print(lines)
        with open(prefix+suffix+'_agg.csv', 'w') as file1:
            for line in lines: file1.write(line)
        spliced = pd.read_csv(prefix+suffix+'_agg.csv', usecols=range(0, 5))
    return spliced

def splice_by_node(prefix, suffix, num_nodes, extra_overlap=0):
    """
    Use splice_by_max_column for each node, with prefix=prefix+node.
    Then, drop 1+extra_overlap
    Combine all files with matching prefix and suffix, each with the same number of states,
    by assigning rows based on the maximum in the given column.

    :param int num_nodes: The number of nodes to splice.

    >>> import pandas as pd
    >>> from pyfeasst import multistate_accumulator
    >>> from pyfeasst import macrostate_distribution
    >>> spliced = splice_by_node(prefix='../../tests/lj_enn', suffix='.txt', num_nodes=2)
    >>> round(float(spliced['average'][375]), 8)
    -2001.76687973
    >>> round(float(spliced['average'][376]), 8)
    -2012.46764871
    >>> len(spliced)
    476
    >>> spliced.to_csv('spliced.csv')
    >>> lnpi = macrostate_distribution.splice_files(prefix='../../tests/lj_lnpin', suffix='.txt')
    >>> lnpi.concat_dataframe(spliced, add_prefix='e_')
    >>> round(float(lnpi.equilibrium()), 8)
    -0.31402411
    >>> vapor, liquid = lnpi.split()
    >>> round(float(-vapor.ln_prob()[0]*0.7/8**3), 8)  # pressure
    0.00136904
    >>> round(float(vapor.ensemble_average('e_average')/vapor.average_macrostate()), 8)
    -0.02500369
    >>> round(float(liquid.ensemble_average('e_average')/liquid.average_macrostate()), 8)
    -6.09838831
    """
    node_data = list()
    for node in range(num_nodes):
        dat = splice_by_max_column(prefix=prefix+str(node), suffix=suffix)
        if node > 0:
            dat.drop(dat.head(1 + extra_overlap).index, inplace=True)
            dat['state'] += node_data[-1]['state'].values[-1]
        node_data.append(dat)
    data = pd.concat(node_data)
    data.reset_index(inplace=True)
    block = 0
    while 'block'+str(block) in data.columns:
        cname = 'block'+str(block)
        # if data[cname].isnull().values.any():
        #     data.drop([cname], inplace=True, axis=1)
        block += 1
    return data

if __name__ == "__main__":
    import doctest
    doctest.testmod()
