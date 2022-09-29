"""
This module computes coarse-grained quantities from PDB files
"""

import numpy as np
import pandas as pd
from biopandas.pdb import PandasPdb

ELEMENT_MASS = {'H': 1.00784, 'C': 12.011, 'N': 14.0067, 'O': 15.999, 'S': 32.065}
ELEMENT_DIAMETER = {'H': 0.4, 'C': 3.7, 'N': 3.3, 'O': 3.1, 'S': 3.9}

def subset(pdb_file, chains):
    """
    Return the subset of ATOM in PDB file that corresponds to the given residues of each chain.

    :param str pdf_file: the file name of the pdb
    :param dict chains: a dictionary of a chain, followed by the residues in that chain

    >>> from pyfeasst import coarse_grain_pdb
    >>> fc = coarse_grain_pdb.subset(pdb_file="../../tests/1igt.pdb", chains={'B': range(236, 244), 'D': range(236, 244)})
    >>> fc['chain_id'].values[0]
    'B'
    >>> fc['residue_number'].values[0]
    236
    >>> fc['chain_id'].values[-1]
    'D'
    >>> fc['residue_number'].values[-1]
    243
    """
    pdb = PandasPdb().read_pdb(pdb_file)
    subset = list()
    for chain in chains:
        sub = pdb.df['ATOM'].loc[(pdb.df['ATOM']['chain_id'] == chain)]
        sub = sub.loc[(sub['residue_number'].isin(chains[chain]))]
        subset.append(sub)
    return pd.concat(subset)

def center_of_mass(subset):
    """
    Return the center of mass position of the given subset.

    >>> from pyfeasst import coarse_grain_pdb
    >>> hinge = coarse_grain_pdb.subset(pdb_file="../../tests/1igt.pdb", chains={'B': range(236, 244), 'D': range(236, 244)})
    >>> r_com_hinge = coarse_grain_pdb.center_of_mass(hinge)
    >>> fc = coarse_grain_pdb.subset(pdb_file="../../tests/1igt.pdb", chains={'B': range(248, 475), 'D': range(248, 475)})
    >>> r_com_fc = coarse_grain_pdb.center_of_mass(fc)
    >>> fc_hinge = r_com_fc - r_com_hinge
    >>> round(np.sqrt(np.dot(fc_hinge, fc_hinge)), 8)
    42.35587498
    """
    masses = list()
    for atom in subset['element_symbol']:
        masses.append(ELEMENT_MASS[atom])
    subset['mass'] = masses
    total_mass = np.sum(masses)
    com = list()
    for coord in ['x_coord', 'y_coord', 'z_coord']:
        com.append(np.sum(subset[coord]*subset['mass'])/total_mass)
    return np.array(com)

def pdb_to_fstprt(subset, fstprt_file):
    """
    Write a fstprt_file that is an all-atom hard-sphere representation
    of the protein using atomic diameters given by
    Grunberger, Lai, Blanco and Roberts, J. Phs. Chem. B, 117, 763 (2013).

    Include an extra last site of type 6 which is the center of mass position that can be used for
    the position of the HS reference in B2 mayer sampling calculations.

    >>> from pyfeasst import coarse_grain_pdb
    >>> fc = coarse_grain_pdb.subset('../../tests/1igt.pdb', {'B': range(248, 474), 'D': range(248, 474)})
    >>> coarse_grain_pdb.pdb_to_fstprt(fc, '1igt_fc.fstprt')
    """
    myfile = open(fstprt_file, 'w')
    params = {'num_sites': len(subset) + 1}
    myfile.write("""# LAMMPS-inspired data file
# site types are H, C, N, O, S, in order

{num_sites} sites

6 site types

Site Properties

0 sigma 0.4 cutoff 0.4
1 sigma 3.7 cutoff 3.7
2 sigma 3.3 cutoff 3.3
3 sigma 3.1 cutoff 3.1
4 sigma 3.9 cutoff 3.9
5 sigma 0.0 cutoff 0.0

Sites

0 5 0 0 0
""".format(**params))
    r_com = center_of_mass(subset)
    for index, atom in enumerate(subset['element_symbol']):
        if atom == "H":
            atom_type = 0
        elif atom == "C":
            atom_type = 1
        elif atom == "N":
            atom_type = 2
        elif atom == "O":
            atom_type = 3
        elif atom == "S":
            atom_type = 4
        else:
            print('unrecognized atom:', atom)
            assert(False)
        myfile.write(str(index + 1) + " " + str(atom_type) + " " +
                     str(subset['x_coord'].values[index] - r_com[0]) + " " +
                     str(subset['y_coord'].values[index] - r_com[1]) + " " +
                     str(subset['z_coord'].values[index] - r_com[2]) + "\n")

if __name__ == "__main__":
    import doctest
    doctest.testmod()
