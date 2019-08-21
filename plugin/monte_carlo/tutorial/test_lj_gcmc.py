import sys
import unittest
import feasst
sys.path.insert(0, feasst.install_dir() + '/plugin/monte_carlo/py') # lj directory
import lj

class TestLJ_GCMC(unittest.TestCase):
    def test(self):
        feasst.seed_random_by_date()
        mc = feasst.MonteCarlo()
        mc.set(lj.system(box_length=8))
        mc.set(feasst.MakeCriteriaMetropolis(feasst.args(
          {"beta": str(1./1.5),
           "chemical_potential": "-8.352321"})))
        mc.add(feasst.MakeTrialTranslate(feasst.args(
          {"weight": "1.", "tunable_param": "2."})))
        feasst.add_trial_transfer(mc, feasst.args(
          {"weight": "1.", "particle_type": "0"}))
        lj.add_analysis(mc, int(1e3))
        mc.attempt(int(1e4))
        self.assertTrue(mc.system().configuration().num_particles() > 0)

if __name__ == "__main__":
    unittest.main()
