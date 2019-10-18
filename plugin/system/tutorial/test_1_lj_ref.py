"""Test the FEASST LJ potential with 2 and 30 particle configurations."""

import unittest
import feasst
import lj_system

def two_particle_config(displacement):
    """Return a two particle configuration with the given displacement between centers"""
    config = lj_system.configuration(num=2)
    select = feasst.SelectList()
    select.particle(1, config)
    config.displace_particle(select, displacement)
    return config

class TestMonteCarlo0LJRef(unittest.TestCase):
    """The FEASST implementation of the LJ potential is tested against known
    cases.
    """
    def test_two_particle(self):
        """Test the LJ potential against analytical calculation of two
        particles
        """
        displacement = feasst.Position(feasst.args({"x": "1.2345", "y": "0", "z": "0"}))
        system = lj_system.system(two_particle_config(displacement))

        # compute the energy of the system
        system.energy()

        # compute the expected analytical LJ and LRC energies
        enlj = 4*(pow(displacement.coord(0), -12) - pow(displacement.coord(0), -6))
        rcut = system.configuration().model_params().cutoff().value(0)
        enlrc = (8./3.)*feasst.PI*system.configuration().num_particles()**2/ \
            system.configuration().domain().volume()*((1./3.)*rcut**(-9) - rcut**(-3))

        # Compare the analytical results with the FEASST computed energies.
        # The energies of the individual potentials (e.g., LJ and LRC) are stored as profiles with
        # indices based on the order that the potentials were initialized.
        # Thus, profile index 0 refers to LJ while 1 refers to LRC.
        # In addition, the last computed value of the energy of all potentials is also stored.
        self.assertAlmostEqual(enlj, system.stored_energy_profile()[0], 15)
        self.assertAlmostEqual(enlrc, system.stored_energy_profile()[1], 15)
        self.assertAlmostEqual(enlj + enlrc, system.stored_energy(), 15)

    def test_srsw_ref_config(self):
        """Test the LJ potential against a configuration of 30 particles.
        In particular, the 4th configuration of the LJ SRSW reference:
        https://www.nist.gov/mml/csd/chemical-informatics-research-group/lennard-jones-fluid-reference-calculations
        """
        config = lj_system.configuration()
        feasst.FileXYZ().load("../../system/test/data/lj_sample_config_periodic4.xyz", config)
        self.assertEqual(30, config.num_particles())
        system = lj_system.system(config=config)
        system.energy()
        enlj = -16.790321304625856
        enlrc = -0.5451660014945704
        self.assertAlmostEqual(enlj, system.stored_energy_profile()[0], 15)
        self.assertAlmostEqual(enlrc, system.stored_energy_profile()[1], 15)
        self.assertAlmostEqual(enlj + enlrc, system.energy(), 15)

if __name__ == "__main__":
    unittest.main()
