import sys
import math
import unittest
import argparse
import feasst
sys.path.insert(0, feasst.install_dir() + '/plugin/monte_carlo/py')
import lj
import pyfeasst
sys.path.insert(0, feasst.install_dir() + '/plugin/flat_histogram/py')
import lj_fh

class TestLJ_FH_benchmark(unittest.TestCase):
    def test(self):
        #return
        feasst.seed_random_by_date()

        serial = False
        serial = True
        if serial:
            criteria = feasst.MakeCriteriaFlatHistogram(feasst.args(
                {"beta": str(1./1.5),
                 "chemical_potential": "-2.35231"}
            ))
            criteria.set(feasst.MakeMacrostateGrowthExpanded(
                feasst.Histogram(feasst.args({"width": "0.5", "max": "10"})),
                #feasst.args({"soft_max": "10"})
            ))
            criteria.set(feasst.MakeBiasTransitionMatrix(feasst.args(
                {"min_sweeps": "10", "num_steps_to_update": str(int(1e6))}
            )))
            mc2 = lj_fh.mc(criteria=criteria, forcefield='data.dimer')
            #print(criteria.write())

        else:
          # parallel
          #windows=[[0,2],[2,15],[15,25]]
          windows=feasst.window(0, 10, 8, 2)
          #print(windows)
          #windows=pyfeasst.vector_vector_to_list(feasst.window(0, 25, 3, 2))

          import multiprocessing as mp
          processes = [mp.Process(target=lj_fh.mc, args=(x, windows[x][0], windows[x][1])) for x in range(len(windows))]
          #processes = [mp.Process(target=mc, args=(window[0], window[1])) for window in windows]

          for p in processes:
            p.start()

          for p in processes:
            p.join()

if __name__ == "__main__":
    unittest.main()
