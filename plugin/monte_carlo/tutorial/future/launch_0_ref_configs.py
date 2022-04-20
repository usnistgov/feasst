import math
import sys
import subprocess
import unittest
import pandas as pd

# Run a feasst simulation with the given parameters
def run_fst(params):
    with open("launch.txt", "w") as myfile: myfile.write("""
MonteCarlo
Configuration {config_params}
Potential Model LennardJones
Potential VisitModel LongRangeCorrections
ThermoParams beta 1000000
Metropolis
Log file_name lj.csv max_precision true clear_file true
Run num_trials 1
""".format(**params))
    syscode = subprocess.call("~/feasst/build/bin/fst < launch.txt > launch.log", shell=True, executable='/bin/bash')
    if syscode > 0: sys.exit(1)

class TestMonteCarloLJRef(unittest.TestCase):
    def test_two_particle(self):
        """Test the LJ potential against analytical calculation of two particles"""
        params = {"displacement": 1.2345}
        with open("two.xyz", "w") as myfile: myfile.write(
"""2
-1 8 8 8
0 0 0 0
1 0 0 {displacement}""".format(**params))
        run_fst({"config_params": "particle_type0 /feasst/forcefield/lj.fstprt xyz_file two.xyz"})
        df = pd.read_csv('lj.csv')
        self.assertEqual(2, df['p0'][0])

        # compute the expected analytical LJ and LRC energies
        enlj = 4*(params["displacement"]**(-12) - params["displacement"]**(-6))
        rcut = 3 #mc.system().configuration().model_params().select("cutoff").value(0)
        enlrc = (8./3.)*math.pi*2**2/ 8**3*((1./3.)*rcut**(-9) - rcut**(-3))

        # Compare the analytical results with the FEASST computed energies.
        self.assertAlmostEqual(enlj, df['LennardJones'][0], 15)
        self.assertAlmostEqual(enlrc, df['LongRangeCorrections'][0], 15)
        self.assertAlmostEqual(enlj + enlrc, df['energy'][0], 15)

    def test_srsw_ref_config(self):
        """Test the LJ potential against a configuration of 30 particles.
        In particular, the 4th configuration of the LJ SRSW reference:
        https://www.nist.gov/mml/csd/chemical-informatics-research-group/lennard-jones-fluid-reference-calculations
        """
        run_fst({"config_params": "cubic_box_length 8 particle_type0 /feasst/forcefield/lj.fstprt \
                  xyz_file /feasst/plugin/configuration/test/data/lj_sample_config_periodic4.xyz"})
        df = pd.read_csv('lj.csv')
        self.assertEqual(df['p0'][0], 30)
        enlj = -16.790321304625856
        enlrc = -0.5451660014945704
        self.assertAlmostEqual(enlj, df['LennardJones'][0], 15)
        self.assertAlmostEqual(enlrc, df['LongRangeCorrections'][0], 15)
        self.assertAlmostEqual(enlj + enlrc, df['energy'][0], 15)

    def test_srsw_ref_config_triclinic(self):
        """Test the LJ potential against a configuration of 300 particles in a trinclinic cell.
        In particular, the 3th configuration of the triclinic LJ SRSW reference:
        https://www.nist.gov/mml/csd/chemical-informatics-group/lennard-jones-fluid-reference-calculations-non-cuboid-cell
        """
        run_fst({"config_params": "side_length0 10.0 side_length1 9.84807753012208 side_length2 9.64974312607518 \
            xy 1.7364817766693041 xz 2.5881904510252074 yz 0.42863479791864567 \
            particle_type0 /feasst/forcefield/lj.fstprt \
            xyz_file /feasst/plugin/configuration/test/data/lj_triclinic_sample_config_periodic3.xyz"})
        df = pd.read_csv('lj.csv')
        self.assertEqual(df['p0'][0], 300)
        enlj = -505.78567945268367
        enlrc = -29.37186430697248
        self.assertAlmostEqual(enlj, df['LennardJones'][0], 15)
        self.assertAlmostEqual(enlrc, df['LongRangeCorrections'][0], 15)
        self.assertAlmostEqual(enlj + enlrc, df['energy'][0], 15)

unittest.main(argv=[''], verbosity=2, exit=False)
