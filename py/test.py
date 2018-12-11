import unittest
import feasst
import pyfeasst

class TestCorePosition(unittest.TestCase):
  def test(self):
    pos = feasst.Position()
    pos.set_to_origin_3D()
    assert(pos.coord(0) == 0)

if __name__ == "__main__":
  unittest.main()
