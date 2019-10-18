import unittest
import fh
import lj_fh

class TestFlatHistogramLJ(unittest.TestCase):
    """Test flat histogram grand canonical ensemble Monte Carlo simulations"""
    def test_serial_5max(self):
        """Compare the free energies and potential energies with the NIST SRSW
        https://www.nist.gov/programs-projects/nist-standard-reference-simulation-website
        https://mmlapps.nist.gov/srs/LJ_PURE/eostmmc.htm
        """
        criteria = fh.criteria_flathist(macro_min=0, macro_max=5, tmmc=True)
        monte_carlo = lj_fh.monte_carlo(criteria=criteria)
        lnpi_previous = [
            [-18.707570324988800000, 1],
            [-14.037373358321800000, 0.037092307087640365],
            [-10.050312091655200000, 0.03696447428346385],
            [-6.458920624988570000, 0.037746391500313385],
            [-3.145637424988510000, 0.03809721387875822],
            [-0.045677458321876000, 0.03845757460933292]
        ]
        energy_previous = [
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
            energy_analyzer = monte_carlo.analyze(monte_carlo.num_analyzers() - 1)
            energy_accumulator = energy_analyzer.analyze(macro).accumulator()
            stdev = (energy_previous[macro][1]**2 + energy_accumulator.block_stdev()**2)**(1./2.)
            self.assertAlmostEqual(
                energy_previous[macro][0],
                energy_accumulator.average(),
                #criteria.bias().ln_macro_prob().value(macro),
                delta=4*stdev
            )

if __name__ == "__main__":
    unittest.main()
