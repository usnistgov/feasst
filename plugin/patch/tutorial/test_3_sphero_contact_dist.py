"""
Find the contact distance for two hard spherocylinders with a given orientation.
"""

import argparse
import json
import subprocess
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation as R
from scipy.optimize import minimize
from pyfeasst import fstio
from make_spherocylinder import hard_spherocylinder

# Parse arguments from command line or change their default values.
PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--feasst_install', type=str, default='../../../build/',
                    help='FEASST install directory (e.g., the path to build)')
PARSER.add_argument('--cylinder_length', type=float, default=3, help='large cylinder length (distance between center of end caps)')
PARSER.add_argument('--cylinder_diameter', type=float, default=1, help='large cylinder diameter')
PARSER.add_argument('--cubic_side_length', type=int, default=100, help='cubic periodic boundary length')
PARSER.add_argument('--num_orientations_per_half_pi', type=int, default=2, help='number of orientations per 90 degrees')

# Convert arguments into a parameter dictionary, and add argument-dependent parameters.
ARGS, UNKNOWN_ARGS = PARSER.parse_known_args()
assert len(UNKNOWN_ARGS) == 0, 'An unknown argument was included: '+str(UNKNOWN_ARGS)
PARAMS = vars(ARGS)
PARAMS['script'] = __file__
PARAMS['prefix'] = 'contact'

hard_spherocylinder(length=PARAMS['cylinder_length'],
                    diameter=PARAMS['cylinder_diameter'],
                    file_name=PARAMS['prefix'] + '.fstprt')

# write xyz relative orientation of two solids of revolution as described in feasst::SolidOfRevolutionTable
# place the center of both solids on the origin
def write_xyz(costheta1, costheta2, cospsi):
    i = np.array([0, 0, 1])
    j = np.array([0, 0, 1])
    r = R.from_euler('y', 0.5*np.pi - np.arccos(costheta1))
    i = r.apply(i)
    r = R.from_euler('y', -0.5*np.pi + np.arccos(costheta2))
    j = r.apply(j)
    r = R.from_euler('x', np.arccos(cospsi))
    j = r.apply(j)
#    j += np.array([separation, 0., 0.])
    xyz_params = {'cubic_side_length': PARAMS['cubic_side_length'],
                  'x1': i[0], 'y1': i[1], 'z1': i[2],
                  'x3': j[0], 'y3': j[1], 'z3': j[2]}
    with open(PARAMS['prefix'] + '_init.xyz', 'w') as f:
        f.write("""4
-1 {cubic_side_length} {cubic_side_length} {cubic_side_length} 0 0 0
0 0 0 0
1 {x1} {y1} {z1}
0 0 0 0
1 {x3} {y3} {z3}""".format(**xyz_params))

def write_feasst_script(params, script_file):
    """ Write fst script for a single simulation with keys of params {} enclosed. """
    with open(script_file, 'w', encoding='utf-8') as myfile:
        myfile.write("""
MonteCarlo
Configuration cubic_side_length {cubic_side_length} particle_type0 {prefix}.fstprt \
    group0 centers centers_site_type0 0 \
    add_particles_of_type0 2 xyz_file {prefix}_init.xyz
Potential Model HardSphere VisitModelInner Spherocylinder group centers
TwoParticleContact output_file {prefix}_dist.txt tolerance 1e-6 dimension 0
""".format(**params))

def get_contact(costheta1, costheta2, cospsi):
    write_xyz(costheta1=costheta1, costheta2=costheta2, cospsi=cospsi)
    #syscode = subprocess.call(PARAMS['feasst_install'] + 'bin/fst < ' + PARAMS['prefix'] + '.txt > '+PARAMS['prefix'] + '_run.log', shell=True, executable='/bin/bash')
    syscode = subprocess.call(PARAMS['feasst_install'] + 'bin/fst < ' + PARAMS['prefix'] + '.txt > /dev/null 2>&1', shell=True, executable='/bin/bash')
    df = pd.read_csv(PARAMS['prefix'] + '_dist.txt', header=None)
    #print(sep, df['energy'][0])
    return df[0][0]

write_feasst_script(PARAMS, PARAMS['prefix']+'.txt')

#def objective_function(sep):
#    dist = sep[0]
#    #print('dist', dist)
#    if dist < 0.5*PARAMS['cylinder_diameter']:
#      return 1000
#    en = get_en(dist, costheta1, costheta2, cospsi)
#    #print('dist', dist, 'en', en)
#    if en > 1:
#      return 10000-dist
#    if en == 0:
#      return 100+dist

dcs = 1./(PARAMS['num_orientations_per_half_pi'])
contacts = list()
#print(get_contact(0, 0, 1))
for costheta1 in np.arange(0, 1+1e-8, dcs):
    for costheta2 in np.arange(0, 1+1e-8, dcs):
        for cospsi in np.arange(-1, 1+1e-8, dcs):
            contacts.append(get_contact(costheta1, costheta2, cospsi))
            #res = minimize(objective_function, [1.5], method='Nelder-Mead', tol=1e-6)
            #print(res)
            #contacts.append(res.x[0])
#print(contacts)

table_file = 'tablek'+str(PARAMS['num_orientations_per_half_pi'])+'l'+str(PARAMS['cylinder_length'])+'d'+str(PARAMS['cylinder_diameter'])+'.txt'
with open(table_file, 'w') as file1:
    file1.write("site_types 1 1\nnum_orientations_per_half_pi " + str(PARAMS['num_orientations_per_half_pi'])+"\ngamma 1\ndelta 0\nnum_z 0\nsmoothing_distance 0\n")
    for _, con in enumerate(contacts):
        file1.write(str(round(con, 6))+'\n')
