"""
This module computes coarse-grained quantities from PDB files
"""

import numpy as np
import pandas as pd
from biopandas.pdb import PandasPdb
pd.options.mode.chained_assignment = None

ELEMENT_MASS = {'H': 1.00784, 'C': 12.011, 'N': 14.0067, 'O': 15.999, 'S': 32.065}
ELEMENT_DIAMETER = {'H': 0.4, 'C': 3.7, 'N': 3.3, 'O': 3.1, 'S': 3.9}

def subset(pdb_file, chains, skip_hydrogens=False):
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
        if skip_hydrogens:
            sub = pdb.df['ATOM'].loc[(pdb.df['ATOM']['chain_id'] == chain) &
                                     (pdb.df['ATOM']['element_symbol'] != 'H')]
        else:
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

def pdb_to_fstprt(subset, fstprt_file, add_com=True, skip_hydrogens=False):
    """
    Write a fstprt_file that is an all-atom hard-sphere representation
    of the protein using atomic diameters given by
    Grunberger, Lai, Blanco and Roberts, J. Phs. Chem. B, 117, 763 (2013).

    If add_com, include an extra site which is the center of mass position that can be used for
    the position of the HS reference in B2 mayer sampling calculations.

    >>> from pyfeasst import coarse_grain_pdb
    >>> fc = coarse_grain_pdb.subset('../../tests/1igt.pdb', {'B': range(248, 474), 'D': range(248, 474)})
    >>> coarse_grain_pdb.pdb_to_fstprt(fc, '1igt_fc.fstprt')
    """
    myfile = open(fstprt_file, 'w')
    params = {'num_sites': len(subset) + 1}
    params['htype'] = '4 sigma 0.4 cutoff 0.4'
    params['num_site_types'] = 5
    if skip_hydrogens:
        params['htype'] = ''
        params['num_site_types'] -= 1
    if add_com:
        params['num_site_types'] += 1
        if skip_hydrogens:
            params['htype'] = '4 sigma 0.0 cutoff 0.0'
            params['com_prop'] = ''
            params['com_type'] = '4'
        else:
            params['com_prop'] = '5 sigma 0.0 cutoff 0.0'
            params['com_type'] = '5'
    myfile.write("""# LAMMPS-inspired data file
# site types are H, C, N, O, S, in order

{num_sites} sites

{num_site_types} site types

Site Properties

0 sigma 3.7 cutoff 3.7
1 sigma 3.3 cutoff 3.3
2 sigma 3.1 cutoff 3.1
3 sigma 3.9 cutoff 3.9
{htype}
{com_prop}

Sites

0 {com_type} 0 0 0
""".format(**params))
    r_com = center_of_mass(subset)
    for index, atom in enumerate(subset['element_symbol']):
        if atom == "H":
            atom_type = 4
        elif atom == "C":
            atom_type = 0
        elif atom == "N":
            atom_type = 1
        elif atom == "O":
            atom_type = 2
        elif atom == "S":
            atom_type = 3
        else:
            print('unrecognized atom:', atom)
            assert(False)
        myfile.write(str(index + 1) + " " + str(atom_type) + " " +
                     str(subset['x_coord'].values[index] - r_com[0]) + " " +
                     str(subset['y_coord'].values[index] - r_com[1]) + " " +
                     str(subset['z_coord'].values[index] - r_com[2]) + "\n")

def write_fstprt(mol):
    """
    Write a fstprt file that is an all-atom representation of the protein using
    AutoDock parameters.

    :param str mol: the name of the mol.pqr file.
    """
    myfile = open(mol+".fstprt.log", 'w')
    # list the epsilon (kcal/mol) and radius, r (A) for each atom type, and alternate labels
    atoms = {
      'H': {'eps': 0.02, 'rad': 0.8909, 'mass': 1.00784, 'alt': ['H', 'H2', 'H3', 'HH', 'HH1', 'HH2', 'HA', 'HA1', 'HA2', 'HA3', 'HB', 'HB1', 'HB2', 'HB3', 'HG', 'HG1', 'HG2', 'HG11', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23', 'HG2', 'HG3', 'HD1', 'HD11', 'HD12', 'HD13', 'HD21', 'HD22', 'HD23', 'HD2', 'HD3', 'HE', 'HE1', 'HE2', 'HE21', 'HE22', 'HE23', 'HE3', 'HH11', 'HH12', 'HH21', 'HH22', 'HZ', 'HZ1', 'HZ2', 'HZ3']},
      'C': {'eps': 0.15, 'rad': 1.7818, 'mass': 12.011, 'alt': ['C', 'CA', 'CB', 'CG', 'CG1', 'CG2', 'CD', 'CD1', 'CD2', 'CE', 'CE1', 'CE2', 'CE3', 'CZ', 'CZ1', 'CZ2', 'CZ3', 'CH', 'CH1', 'CH2']},
      'N': {'eps': 0.16, 'rad': 1.5591, 'mass': 14.0067, 'alt': ['N', 'NH1', 'NH2', 'ND1', 'ND2', 'NE', 'NE1', 'NE2', 'NZ']},
      'O': {'eps': 0.20, 'rad': 1.4254, 'mass': 15.999, 'alt': ['O', 'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1', 'OG2', 'OH', 'OXT']},
      'S': {'eps': 0.20, 'rad': 1.7818, 'mass': 32.065, 'alt': ['S', 'SG', 'SD']},
    }

    df = pd.read_csv(mol+'.pqr', header=None, sep='\s+')
    df = df[df[0] != 'TER']
    df = df[:len(df)-1]
    df.reset_index(inplace=True)
    df['atom'] = df[2]
    df['mass'] = df[1]*0
    df['sig'] = df[1]*0
    df['eps'] = df[1]*0

    #print(df[2])
    #print(df[2].values)

    for _, atom in enumerate(atoms):
        #print(atom)
        #print(atoms[atom]['eps'])
        for _, aa in enumerate(atoms[atom]['alt']):
            df.loc[df[2]==aa, 'atom'] = atom
            #df['atom'][df[2]==aa] = atom
            df.loc[df[2]==aa, 'mass'] = atoms[atom]['mass']
            #df['mass'][df[2]==aa] = atoms[atom]['mass']
            df.loc[df[2]==aa, 'eps'] = atoms[atom]['eps']
            #df['eps'][df[2]==aa] = atoms[atom]['eps']
            df.loc[df[2]==aa, 'sig'] = 2*atoms[atom]['rad']
            #df['sig'][df[2]==aa] = 2*atoms[atom]['rad']

    assert len(df[df['sig'] == 0]) == 0
    #print(df[df['sig'] == 0])
    #print(df['eps'])
    df.to_csv(mol+'.vdw')

    # net mass
    mass = np.sum(df['mass'])
    myfile.write('mass ' + str(mass))

    # v_scale
    m_kda = mass/1e3
    v_s = 0.18*(0.27*m_kda+80)/(m_kda + 80)
    myfile.write('v_s ' + str(v_s))
    v_s *= 4.184 # convert kcal/mol to kJ/mol

    # center of mass
    cmx = np.sum(df[5]*df['mass'])/np.sum(df['mass'])
    cmy = np.sum(df[6]*df['mass'])/np.sum(df['mass'])
    cmz = np.sum(df[7]*df['mass'])/np.sum(df['mass'])
    myfile.write('cmx ' + str(cmx))
    myfile.write('cmy ' + str(cmy))
    myfile.write('cmz ' + str(cmz))

    # make a unique type for each possible combination of charge, epsilon and sigma
    atom_type = list()
    type_index = 0
    df['type'] = 0
    for atom in range(len(df)):
        #print('atom', atom)
        q = df[8][atom]
        e = df['eps'][atom]
        s = df['sig'][atom]
        new = True
        for iat, at in enumerate(atom_type):
            #print('at', at)
            if at['q'] == q and at['sig'] == s and at['eps'] == e:
                new = False
                df.loc[atom, 'type'] = iat
                #df['type'][atom] = iat
        if new:
            df.loc[atom, 'type'] = len(atom_type)
            #df['type'][atom] = len(atom_type)
            atom_type.append({'q': q, 'sig': s, 'eps': e})

    #print(atom_type)
    #print(df)

    df.to_csv(mol+'_typed.vdw')
    #print(df[5])

    # write a fstprt file
    myfile = open(mol+'.fstprt', 'w')
    params = {
      'num_sites': len(df),
      'num_site_types': len(atom_type),
    }
    myfile.write("""# FEASST particle file

{num_sites} sites

{num_site_types} site types

Site Properties

""".format(**params))
    for iat, at in enumerate(atom_type):
        myfile.write(str(iat)+ " charge "  + str(at['q']) + " epsilon " + str(at['eps']*v_s) + " sigma " + str(at['sig']) + ' cutoff 200\n')

    myfile.write("\n\nSites\n\n")

    for atom in range(len(df)):
        myfile.write(str(atom) + " " + str(df['type'][atom]) + " " + str(df[5][atom]-cmx) + " " + str(df[6][atom]-cmy) + " " + str(df[7][atom]-cmz) + '\n')

    # write a fstprt file with a reference site at the beginning for mayer sampling
    myfile = open(mol+'_with_ref.fstprt', 'w')
    params = {
      'num_sites': len(df)+1,
      'num_site_types': len(atom_type)+1,
    }
    myfile.write("""# FEASST particle file

{num_sites} sites

{num_site_types} site types

Site Properties

0 charge 0 epsilon 0 sigma 0 cutoff 0
""".format(**params))
    for iat, at in enumerate(atom_type):
        myfile.write(str(iat+1)+ " charge "  + str(at['q']) + " epsilon " + str(at['eps']*v_s) + " sigma " + str(at['sig']) + ' cutoff 200\n')

    myfile.write("\n\nSites\n\n")

    myfile.write("0 0 0 0 0\n")
    for atom in range(len(df)):
        myfile.write(str(atom+1) + " " + str(df['type'][atom]+1) + " " + str(df[5][atom]-cmx) + " " + str(df[6][atom]-cmy) + " " + str(df[7][atom]-cmz) + '\n')
if __name__ == "__main__":
    import doctest
    doctest.testmod()
