import sys
import math
import unittest
import argparse
import feasst
sys.path.insert(0, feasst.install_dir() + '/plugin/monte_carlo/py') # lj directory
import lj
import pyfeasst
import test_lj_fh

class TestLJ_FH_benchmark(unittest.TestCase):
    def test(self):
        #return
        feasst.seed_random_by_date()

        #windows=[[0,2],[2,15],[15,25]]
        windows=feasst.window(0, 370, 8, 2)
        #print(windows)
        #windows=pyfeasst.vector_vector_to_list(feasst.window(0, 25, 3, 2))

        serial = False
        serial = True
        if serial:
            criteria, mc2 = test_lj_fh.mc(proc=0, macro_min=0, macro_max=370, tmmc=True)
            #print(criteria.write())

        else:
          # parallel
          import multiprocessing as mp
          processes = [mp.Process(target=test_lj_fh.mc, args=(x, windows[x][0], windows[x][1])) for x in range(len(windows))]
          #processes = [mp.Process(target=mc, args=(window[0], window[1])) for window in windows]

          for p in processes:
            p.start()

          for p in processes:
            p.join()

if __name__ == "__main__":
    unittest.main()
