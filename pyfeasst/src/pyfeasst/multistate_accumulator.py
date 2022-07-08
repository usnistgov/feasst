"""
This module analyzes and manipulates multistate accumulators.
For example, the canonical ensemble average energy for a given number of particles.
"""

from pathlib import Path
import pandas as pd

def splice_by_max_column(prefix, suffix, column='moment0'):
    """
    Combine all files with matching prefix and suffix, each with the same number of states,
    by assigning rows based on the maximum in the given column.

    :param str prefix: Find all files beginning with this prefix.
    :param str suffix: Find all files endding with this suffix.
    :param str column:
        Use the row with the maximum column. For example, moment0 is the highest number of samples.

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
    for filename in Path('.').rglob(prefix+'*'+suffix):
        if spliced is None:
            spliced = pd.read_csv(filename)
        else:
            addition = pd.read_csv(filename)
            newer = addition[addition[column] > spliced[column]]
            if len(newer) > 0:
                spliced.loc[newer['state'].values[0]: newer['state'].values[-1]] = newer
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
    >>> spliced['average'][375]
    -2001.7668797287809
    >>> spliced['average'][376]
    -2012.4676487137983
    >>> len(spliced)
    476
    >>> spliced.to_csv('spliced.csv')
    >>> lnpi = macrostate_distribution.splice_files(prefix='../../tests/lj_lnpin', suffix='.txt')
    >>> lnpi.concat_dataframe(spliced, add_prefix='e_')
    >>> lnpi.equilibrium()
    -0.31402410845849854
    >>> vapor, liquid = lnpi.split()
    >>> -vapor.ln_prob()[0]*0.7/8**3  # pressure
    0.0013690413500822545
    >>> vapor.ensemble_average('e_average')/vapor.average_macrostate()
    -0.025003690574197115
    >>> liquid.ensemble_average('e_average')/liquid.average_macrostate()
    -6.098388307322735
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
