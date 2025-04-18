"""
Generate a PSF file for a number of united atom linear chains.
"""

import argparse

PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
PARSER.add_argument('--num_particles', type=int, default=1, help='number of particles')
PARSER.add_argument('--sites_per_particle', type=int, default=8, help='number of sites in a linear particle')
PARSER.add_argument('--output_psf', '-o', type=str, default='n-octane.psf', help='output file name')
ARGS, UNKNOWN_ARGS = PARSER.parse_known_args()
assert len(UNKNOWN_ARGS) == 0, 'An unknown argument was included: '+str(UNKNOWN_ARGS)
PARAMS = vars(ARGS)

PARAMS['num_sites'] = PARAMS['num_particles']*PARAMS['sites_per_particle']
psf = open(PARAMS['output_psf'], 'w')
psf.write("""PSF CMAP

       {num_sites} !NATOM\n""".format(**PARAMS))
PARAMS['site'] = 1
for part in range(PARAMS['num_particles']):
    for site in range(PARAMS['sites_per_particle']):
        if site == 0 or site == PARAMS['sites_per_particle'] - 1:
            PARAMS['element'] = 'N'
        else:
            PARAMS['element'] = 'C'
        psf.write("""{site:8d} U    1    MET  {element}    NH3   -0.300000       14.0070           0\n""".format(**PARAMS))
        PARAMS['site'] += 1

PARAMS['num_bonds'] = PARAMS['num_particles']*(PARAMS['sites_per_particle'] - 1)
psf.write("""
       {num_bonds} !NBOND: bonds""".format(**PARAMS))

PARAMS['bond'] = 1
PARAMS['site'] = 1
for part in range(PARAMS['num_particles']):
    for site in range(PARAMS['sites_per_particle']):
        if site < PARAMS['sites_per_particle'] - 1:
            PARAMS['next_site'] = PARAMS['site'] + 1
            if PARAMS['bond'] % 4 == 1:
                psf.write("\n")
            psf.write("""{site:8d}{next_site:8d}""".format(**PARAMS))
            PARAMS['bond'] += 1
        PARAMS['site'] += 1
psf.close()
