import sys
import unittest
import feasst
import pyfeasst
sys.path.insert(0, feasst.install_dir() + '/plugin/flat_histogram/py')
import lj_fh

class TestLJ_FH_benchmark(unittest.TestCase):
    def test(self):
        feasst.seed_random_by_date()
        # serial = False
        serial = True
        if serial:
            criteria = lj_fh.criteria_flathist(macro_min=0, macro_max=370, tmmc=True)
            mc2 = lj_fh.mc(criteria=criteria)
            #print(criteria.write())

        else:
            # parallel
            #windows=[[0,2],[2,15],[15,25]]
            num_procs = 3
            windows=feasst.window(0, 370, num_procs, 2)
            print(windows)
            #windows=pyfeasst.vector_vector_to_list(feasst.window(0, 25, 3, 2))

            # HWH issue with multiprocessing if criterium object is provided as arg
            import multiprocessing as mp
  #          with mp.Pool(processes = num_procs) as pool:
  #              pool.starmap(lj_fh.mc, zip(range(num_procs), repeat(criteria=lj_fh.criteria_flathist(macro_min=

            processes = []
          # criterium = []
            for x in range(len(windows)):
              # criterium.append(lj_fh.criteria_flathist(macro_min=windows[x][0], macro_max=windows[x][1]))
              # processes.append(mp.Process(target=lj_fh.mc, args=(x,)))
                processes.append(mp.Process(target=lj_fh.mc, args=(x, lj_fh.criteria_flathist(macro_min=windows[x][0], macro_max=windows[x][1]))))
              # processes.append(mp.Process(target=lj_fh.mc, args=(x, criterium[x], )))
            #processes = [mp.Process(target=lj_fh.mc, args=(x,)) for x in range(len(windows))]
            #processes = [mp.Process(target=lj_fh.mc, args=(x, lj_fh.criteria_flathist(windows[x][0], windows[x][1]))) for x in range(len(windows))]
           #processes = [mp.Process(target=lj_fh.mc, args=(x, lj_fh.criteria_flathist(windows[x][0], windows[x][1]))) for x in range(len(windows))]

            for p in processes:
                p.start()

            for p in processes:
                p.join()

if __name__ == "__main__":
    unittest.main()
