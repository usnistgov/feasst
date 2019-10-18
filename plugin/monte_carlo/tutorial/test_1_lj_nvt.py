import sys
import unittest
import feasst
import analyze
sys.path.insert(0, feasst.install_dir() + '/plugin/system/tutorial/') # lj_system
import lj_system

class TestMonteCarlo1LJNVT(unittest.TestCase):
    """Test a canonical ensemble Lennard Jones Monte Carlo simulation"""
    def test_srsw(self, num_particles=500, density=0.001, steps_per=1e5):
        """Compare with the reported average energy from the NIST SRSW.
        https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm
        https://www.nist.gov/programs-projects/nist-standard-reference-simulation-website

        num_particles -- number of LJ particles
        density -- number density
        steps_per -- steps between each Anaylze/Modify
        """
        monte_carlo = feasst.MonteCarlo()
        monte_carlo.set(lj_system.system(lj_system.configuration(
            box_length=(num_particles/density)**(1./3.))))
        monte_carlo.set(feasst.MakeCriteriaMetropolis(feasst.args(
            {"beta": str(1./0.9),
             "chemical_potential": "1."})))
        monte_carlo.add(feasst.MakeTrialTranslate(feasst.args(
            {"weight": "1.", "tunable_param": "2."})))
        monte_carlo.seek_num_particles(num_particles)
        analyze.add(monte_carlo, steps_per)

        # equilibrate
        monte_carlo.attempt(int(1e7))

        # compute average energy using a stepper/analysis and output into file
        energy = feasst.MakeEnergy(feasst.args({
            "steps_per_update": "1",
            "steps_per_write": str(steps_per),
            "file_name": "lj_nvt_srsw_energy.txt",
        }))
        monte_carlo.add(energy)

        # compute average using this script
        energy_alt = feasst.Accumulator()

        # production
        for _ in range(int(1e7)):
            monte_carlo.attempt(1)
            energy_alt.accumulate(monte_carlo.criteria().current_energy())

        # test the average against the NIST SRSW
        stdev = (energy.energy().block_stdev()**2 + (1.89E-05)**2)**(1./2.)
        self.assertAlmostEqual(-9.9165E-03*num_particles, energy.energy().average(),
                               delta=2.576*stdev)

        # test that the two methods of compute average energy are the same
        self.assertAlmostEqual(energy.energy().average(), energy_alt.average(), delta=1e-6)

if __name__ == "__main__":
    unittest.main()
