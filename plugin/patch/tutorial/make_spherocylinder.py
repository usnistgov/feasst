# write particle file
def hard_spherocylinder(length, diameter, file_name):
    with open(file_name, 'w') as f:
        fstprt_params = {'central_cut': length + diameter, 'length': length, 'diameter': diameter}
        f.write("""# LAMMPS-inspired data file

2 sites
1 bonds

2 site types
1 bond types

Site Properties

0 sigma {diameter} epsilon 1.0 cutoff {central_cut}
1 sigma {diameter} epsilon 1.0 cutoff {diameter} director 1 spherocylinder_length {length}

Bond Properties

0 RigidBond length 1.0 delta 0.000001

Sites

0 0 0.0 0.0 0.0
1 1 1.0 0.0 0.0

Bonds

0 0 0 1""".format(**fstprt_params))

