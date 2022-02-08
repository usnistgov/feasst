import sys
import subprocess
import unittest
import feasst as fst
import pyfeasst

# Run a feasst simulation with the given parameters
def run_fst(params):
    with open("launch.txt", "w") as myfile: myfile.write("""
Configuration {config_params}
Potential Model LennardJones
Potential VisitModel LongRangeCorrections
ThermoParams
Metropolis
Checkpoint file_name checkpoint.fst
WriteCheckpoint
""".format(**params))
    syscode = subprocess.call(fst.install_dir() + "/build/bin/fst < launch.txt > launch.log", shell=True, executable='/bin/bash')
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
        mc = fst.MonteCarlo().deserialize(pyfeasst.read_checkpoint("checkpoint.fst"))
        self.assertEqual(mc.configuration().num_particles(), 2)

        # compute the expected analytical LJ and LRC energies
        enlj = 8*(params["displacement"]**(-12) - params["displacement"]**(-6))
        rcut = mc.system().configuration().model_params().select("cutoff").value(0)
        enlrc = (8./3.)*fst.PI*mc.system().configuration().num_particles()**2/ \
            mc.configuration().domain().volume()*((1./3.)*rcut**(-9) - rcut**(-3))

        # Compare the analytical results with the FEASST computed energies.
        # The energies of the individual potentials (e.g., LJ and LRC) are stored as profiles with
        # indices based on the order that the potentials were initialized.
        # Thus, profile index 0 refers to LJ while 1 refers to LRC.
        # In addition, the last computed value of the energy of all potentials is also stored.
        self.assertAlmostEqual(enlj, mc.system().stored_energy_profile()[0], 15)
        self.assertAlmostEqual(enlrc, mc.system().stored_energy_profile()[1], 15)
        self.assertAlmostEqual(enlj + enlrc, mc.system().stored_energy(), 15)

    def test_srsw_ref_config(self):
        """Test the LJ potential against a configuration of 30 particles.
        In particular, the 4th configuration of the LJ SRSW reference:
        https://www.nist.gov/mml/csd/chemical-informatics-research-group/lennard-jones-fluid-reference-calculations
        """
        run_fst({"config_params": "cubic_box_length 8 particle_type0 /feasst/forcefield/lj.fstprt \
                  xyz_file /feasst/plugin/configuration/test/data/lj_sample_config_periodic4.xyz"})
        mc = fst.MonteCarlo().deserialize(pyfeasst.read_checkpoint("checkpoint.fst"))
        self.assertEqual(mc.configuration().num_particles(), 30)
        enlj = -16.790321304625856
        enlrc = -0.5451660014945704
        self.assertEqual(enlj, mc.system().stored_energy_profile()[0], 15)
        self.assertAlmostEqual(enlrc, mc.system().stored_energy_profile()[1], 15)
        self.assertAlmostEqual(enlj + enlrc, mc.system().energy(), 15)

    def test_srsw_ref_config_triclinic(self):
        """Test the LJ potential against a configuration of 300 particles in a trinclinic cell.
        In particular, the 3th configuration of the triclinic LJ SRSW reference:
        https://www.nist.gov/mml/csd/chemical-informatics-group/lennard-jones-fluid-reference-calculations-non-cuboid-cell
        """
        run_fst({"config_params": "side_length0 10.0 side_length1 9.84807753012208 side_length2 9.64974312607518 \
            xy 1.7364817766693041 xz 2.5881904510252074 yz 0.42863479791864567 \
            particle_type0 /feasst/forcefield/lj.fstprt \
            xyz_file /feasst/plugin/configuration/test/data/lj_triclinic_sample_config_periodic3.xyz"})
        mc = fst.MonteCarlo().deserialize(pyfeasst.read_checkpoint("checkpoint.fst"))
        self.assertEqual(mc.configuration().num_particles(), 300)
        enlj = -505.78567945268367
        enlrc = -29.37186430697248
        self.assertEqual(enlj, mc.system().stored_energy_profile()[0], 15)
        self.assertAlmostEqual(enlrc, mc.system().stored_energy_profile()[1], 15)
        self.assertAlmostEqual(enlj + enlrc, mc.system().energy(), 15)

unittest.main(argv=[''], verbosity=2, exit=False)
