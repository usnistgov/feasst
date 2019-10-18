import sys
import unittest
import feasst
import analyze
sys.path.insert(0, feasst.install_dir() + '/plugin/system/tutorial/') # lj_system
import lj_system

class TestMonteCarlo2LJGCMC(unittest.TestCase):
    """Test a grand canonical ensemble Lennard Jones Monte Carlo simulation"""
    def test(self):
        """Compute the average number of particles and assert that it is greater than 0"""
        monte_carlo = feasst.MonteCarlo()
        monte_carlo.set(lj_system.system(lj_system.configuration(box_length=8)))
        monte_carlo.set(feasst.MakeCriteriaMetropolis(feasst.args(
            {"beta": str(1./1.5), "chemical_potential": "-8.352321"})))
        monte_carlo.add(feasst.MakeTrialTranslate(feasst.args(
            {"weight": "1.", "tunable_param": "2."})))
        # add an insertion/deletion trial attempt
        feasst.add_trial_transfer(monte_carlo, feasst.args({"weight": "1.", "particle_type": "0"}))
        steps_per = int(1e3)
        analyze.add(monte_carlo, steps_per)

        # Add an Analyze which computes the average number of particles.
        # Just before adding, store the number of existing Analyzers in order to remember the
        # index of the newly added Analyze.
        analyze_index = monte_carlo.num_analyzers()
        monte_carlo.add(feasst.MakeNumParticles(feasst.args(
            {"steps_per_write": str(steps_per), "file_name": "gcmc_num_particles.txt"})))

        # peform a short simulation
        monte_carlo.attempt(int(1e5))

        # assert that particles were added during the simulation
        self.assertTrue(monte_carlo.analyze(analyze_index).accumulator().average() > 0)

if __name__ == "__main__":
    unittest.main()
