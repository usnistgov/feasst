import sys
import subprocess
import unittest
import pandas as pd

class TestMonteCarlo1LJNVT(unittest.TestCase):
    """Test a canonical ensemble Lennard Jones Monte Carlo simulation"""
    def test_srsw_alt(self, num_particles=500, density=0.001, trials_per=1e5, beta=1./0.9,
                      fstprt="/feasst/forcefield/lj.fstprt",
                      num_equil=1e7, num_prod=1e7):
        """Compare with the reported average energy from the NIST SRSW.
        https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm
        https://www.nist.gov/programs-projects/nist-standard-reference-simulation-website

        num_particles -- number of LJ particles
        density -- number density
        trials_per -- steps between each Anaylze/Modify
        """
        params = {"box_length": (num_particles/density)**(1./3.)}
        params = dict(locals(), **params)
        with open('launch_1_lj_nvt.txt', 'w') as fsttxt:
            fsttxt.write("""
MonteCarlo
Checkpoint file_name checkpoint.fst
RandomMT19937 seed time
Configuration cubic_box_length {box_length} particle_type {fstprt}
Potential Model LennardJones
Potential VisitModel LongRangeCorrections
ThermoParams beta 0.1 chemical_potential 10
Metropolis
TrialTranslate tunable_param 2. tunable_target_acceptance 0.2
TrialAdd particle_type 0
Run until_num_particles {num_particles}
RemoveTrial name TrialAdd
ThermoParams beta {beta}
Tune trials_per {trials_per}
CheckEnergy trials_per {trials_per} tolerance 1e-8

# equilibrate
Run num_trials {num_equil}
RemoveModify name Tune

# production analysis and output
Log trials_per {trials_per} file_name lj.txt
Energy trials_per_write {trials_per} file_name en.txt
Run num_trials {num_prod}
""".format(**params))
        syscode = subprocess.call("~/feasst/build/bin/fst < launch_1_lj_nvt.txt >> launch.log", shell=True, executable='/bin/bash')
        if syscode > 0: sys.exit(1)

        # test the average energy against the NIST SRSW
        df = pd.read_csv('en.txt')
        stdev = (df['block_stdev'][0]**2 + (1.89E-05)**2)**(1./2.)
        self.assertAlmostEqual(-9.9165E-03*num_particles,
                               df['average'][0],
                               delta=2.576*stdev)

unittest.main(argv=[''], verbosity=2, exit=False)
