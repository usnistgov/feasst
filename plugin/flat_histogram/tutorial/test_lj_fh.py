import sys
import math
import unittest
import argparse
import feasst
#sys.path.insert(0, feasst.install_dir() + '/plugin/monte_carlo/py')
#import lj
sys.path.insert(0, feasst.install_dir() + '/plugin/flat_histogram/py')
import lj_fh
import pyfeasst

class TestLJ_WL(unittest.TestCase):
    def test_serial_5max(self):
        feasst.seed_random_by_date()
        criteria = lj_fh.criteria_flathist(macro_min=0, macro_max=5, tmmc=True)
        mc2 = lj_fh.mc()
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
