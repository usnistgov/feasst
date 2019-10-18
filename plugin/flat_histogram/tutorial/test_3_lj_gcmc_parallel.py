import sys
import unittest
import feasst
import pyfeasst
import fh
import lj_fh

class TestFlatHistogramLJGCMCparallel(unittest.TestCase):
    def test(self):
        serial = False
        # serial = True
        if serial:
            criteria = fh.criteria_flathist(macro_min=0, macro_max=370, tmmc=True)
            mc2 = lj_fh.monte_carlo(criteria=criteria)
            #print(criteria.write())

        else:
            # parallel
            #windows=[[0,2],[2,15],[15,25]]
            num_procs = 3
            windows=feasst.window(0, 370, num_procs, 2)
            print(windows)
            #windows=pyfeasst.vector_vector_to_list(feasst.window(0, 25, 3, 2))

            criterium = list()
            for _, winx in enumerate(windows):
                criterium.append(fh.criteria_flathist(macro_min=winx[0], macro_max=winx[1]))

            # HWH issue with multiprocessing if criterium object is provided as arg
            import multiprocessing as mp
#            with mp.Pool(processes = num_procs) as pool:
#                pool.starmap(lj_fh.monte_carlo, zip(range(num_procs), iter(criterium)))
#                repeat(criteria=lj_fh.criteria_flathist(macro_min=

            processes = []
            criterium = list()
            for _, winx in enumerate(windows):
                criterium.append(fh.criteria_flathist(macro_min=winx[0], macro_max=winx[1]))
            for x, crit in enumerate(criterium):
#                #processes.append(mp.Process(target=lj_fh.monte_carlo, args=(x,)))
#                #  processes.append(mp.Process(target=lj_fh.monte_carlo, args=(x, fh.criteria_flathist(macro_min=windows[x][0], macro_max=windows[x][1]))))
                processes.append(mp.Process(target=lj_fh.monte_carlo, args=(x, crit, )))
            #    processes[-1].start()
            #processes = [mp.Process(target=lj_fh.monte_carlo, args=(x,)) for x in range(len(windows))]
            #processes = [mp.Process(target=lj_fh.monte_carlo, args=(x, lj_fh.criteria_flathist(windows[x][0], windows[x][1]))) for x in range(len(windows))]
           #processes = [mp.Process(target=lj_fh.monte_carlo, args=(x, lj_fh.criteria_flathist(windows[x][0], windows[x][1]))) for x in range(len(windows))]

            for p in processes:
                p.start()

            for p in processes:
                p.join()

if __name__ == "__main__":
    unittest.main()
