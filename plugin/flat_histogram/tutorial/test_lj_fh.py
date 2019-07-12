import sys
import math
import unittest
import argparse
import feasst
sys.path.insert(0, feasst.install_dir() + '/plugin/monte_carlo/py') # lj directory
import lj
import pyfeasst

def criteria_flathist(macro_min, macro_max, chemical_potential=-2.352321, tmmc=True):
    criteria = feasst.MakeCriteriaFlatHistogram(feasst.args(
      {"beta": str(1./1.5),
       "chemical_potential": str(chemical_potential)}))
    criteria.set(feasst.MakeMacrostateNumParticles(feasst.Histogram(feasst.args(
      {"width": "1", "min": str(macro_min), "max": str(macro_max)}))))
    if tmmc:
        criteria.set(feasst.MakeBiasTransitionMatrix(feasst.args(
          {"min_sweeps": "20", "num_steps_to_update": str(int(1e6))})))
    else:
        criteria.set(feasst.MakeBiasWangLandau(feasst.args({"min_flatness": "20"})))
    return criteria

def mc(proc, macro_min=0, macro_max=370, steps_per=1e6, tmmc=True):
    mc = feasst.MonteCarlo()
    mc.set(lj.system(box_length=8))

    # add the minimum number of particles to the system
    mc.set(feasst.MakeCriteriaMetropolis(feasst.args(
      {"beta": "0.001", "chemical_potential": "1"})))
    mc.add(feasst.MakeTrialTranslate(feasst.args(
      {"weight": "0.75", "tunable_param": "2."})))
    mc.add(feasst.MakeTrialAdd(feasst.args({"weight": "0.125"})))
    mc.add(feasst.MakeTrialRemove(feasst.args({"weight": "0.125"})))
    mc.seek_num_particles(macro_min)

    # replace the acceptance criteria with flat histogram
    criteria = criteria_flathist(macro_min, macro_max, tmmc=tmmc)
    mc.set(criteria)

    lj.add_analysis(mc, steps_per, proc=proc, log="log"+str(proc)+".txt")
    mc.add(feasst.MakeCriteriaWriter(feasst.args(
      {"steps_per": str(steps_per), "file_name": "crit"+str(proc)+".txt"})))
    mc.add(feasst.MakeEnergy(feasst.args(
      {"file_name": "energy"+str(proc)+".txt",
       "steps_per_update": "1",
       "steps_per_write": str(steps_per),
       "multistate": "true"})))
    mc.add(feasst.MakeCheckpoint(feasst.args(
      {"file_name": "checkpoint"+str(proc)+".txt", "num_hours": "0.01"})))
    mc.run_until_complete()
    #print(mc.criteria().write())
    return criteria, mc

class TestLJ_WL(unittest.TestCase):
    def test_serial_5max(self):
        feasst.seed_random_by_date()
        criteria, mc2 = mc(proc=0, macro_min=0, macro_max=5, tmmc=True)
        lnpi_previous=[
            [-18.707570324988800000, 1],
            [-14.037373358321800000, 0.037092307087640365],
            [-10.050312091655200000, 0.03696447428346385],
            [-6.458920624988570000, 0.037746391500313385],
            [-3.145637424988510000, 0.03809721387875822],
            [-0.045677458321876000, 0.03845757460933292]
        ]
        energy_previous=[
            [0, 1e-14],
            [-0.0006057402333333332, 6.709197666659334e-10],
            [-0.030574223333333334, 9.649146611661053e-06],
            [-0.089928316, 0.0001387472078025413],
            [-0.1784570533333333, 3.3152449884326804e-05],
            [-0.29619201333333334, 1.3487910636322294e-05],
        ]
        for macro in range(criteria.num_states()):
            self.assertAlmostEqual(
                lnpi_previous[macro][0],
                criteria.bias().ln_macro_prob().value(macro),
                delta=lnpi_previous[macro][1]
            )
            energy_accumulator = mc2.analyze(mc2.num_analyzers()-1).analyze(macro).accumulator()
            self.assertAlmostEqual(
                energy_previous[macro][0],
                energy_accumulator.average(),
                #criteria.bias().ln_macro_prob().value(macro),
                delta=2*(energy_previous[macro][1]**2 + energy_accumulator.block_stdev()**2)**(1./2.)
            )

if __name__ == "__main__":
    unittest.main()
