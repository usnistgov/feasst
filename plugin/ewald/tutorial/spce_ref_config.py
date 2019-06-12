import feasst
import unittest
import spce

class TestSPCE_REFCONF(unittest.TestCase):
  def test(self):
    config = feasst.Configuration(feasst.args({"particle_type": "../../../forcefield/data.spce"}))
    feasst.FileXYZ().load("../../system/test/data/spce_sample_config_periodic1.xyz", config)
    system = spce.system(config)
    print(system.energy())
    self.assertAlmostEqual(system.energy(), -4062.47263092246, 10)

if __name__ == "__main__":
    unittest.main()
