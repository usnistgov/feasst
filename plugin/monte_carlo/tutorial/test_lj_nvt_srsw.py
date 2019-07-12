import sys
import unittest
import feasst
sys.path.insert(0, feasst.install_dir() + '/plugin/monte_carlo/py') # lj directory
import lj

class TestLJ_NVT_SRSW(unittest.TestCase):
    def test(self, num_particles=500, density=0.001, steps_per=1e5):
        feasst.seed_random_by_date()
        mc = feasst.MonteCarlo()
        mc.set(lj.system(box_length=(num_particles/density)**(1./3.)))
        mc.set(feasst.MakeCriteriaMetropolis(feasst.args(
          {"beta": str(1./0.9),
           "chemical_potential": "1."})))
        mc.add(feasst.MakeTrialTranslate(feasst.args(
          {"weight": "1.", "tunable_param": "2."})))
        lj.add_analysis(mc, steps_per)
        mc.seek_num_particles(num_particles)

        mc.add(feasst.MakeNumParticles(feasst.args({
          "steps_per_write": str(steps_per),
        })))

        # equilibrate
        mc.attempt(int(1e7))

        # compute average energy using a stepper/analysis
        energy = feasst.MakeEnergy(feasst.args({
          "steps_per_update": "1",
          "steps_per_write": str(steps_per),
          "file_name": "lj_nvt_srsw_energy.txt",
        }))
        mc.add(energy);

        # compute average using this script
        energy_alt = feasst.Accumulator()

        # production
        for trial in range(int(1e7)):
          mc.attempt(1)
          energy_alt.accumulate(mc.criteria().current_energy())

        # test the average against the NIST SRSW
        # https:#mmlapps.nist.gov/srs/LJ_PURE/mc.htm
        # https:#www.nist.gov/programs-projects/nist-standard-reference-simulation-website
        self.assertAlmostEqual(-9.9165E-03*num_particles, energy.energy().average(),
          delta=2.576*(energy.energy().block_stdev()**2 + (1.89E-05)**2)**(1./2.))

        # test that the two averages are the same
        self.assertAlmostEqual(energy.energy().average(), energy_alt.average(), delta=1e-6)

if __name__ == "__main__":
    unittest.main()
