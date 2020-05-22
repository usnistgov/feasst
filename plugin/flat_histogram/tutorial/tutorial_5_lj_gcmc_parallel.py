import sys
import unittest
import multiprocessing as mp
import feasst as fst
import pyfeasst
sys.path.insert(0, fst.install_dir() + '/plugin/flat_histogram/tutorial/')
import fh

class TestFlatHistogramLJGCMCparallel(unittest.TestCase):
    def test(self):
        serial = False
        # serial = True
        if serial:
            criteria = fh.criteria_flathist(macro_min=0, macro_max=370, tmmc=True)
            mc2 = fh.monte_carlo(criteria=criteria)

        else:
            windows=fst.WindowExponential(fst.args({
              "alpha": "2",
              "num": "4",
              "maximum": "50",
              "extra_overlap": "2"})).boundaries()
            #windows=[[0,2],[2,15],[15,25]]
            print(windows)

##            jobs = list()
#            clones = fst.Clones()
#            for proc, win in enumerate(windows):
#                print(proc, win)
#                clones.add(fh.monte_carlo(proc=proc,
#                    criteria=fh.criteria_flathist(macro_min=win[0], macro_max=win[1]),
#                    run=False,
#                    steps_per=1))
#
#            for proc, win in enumerate(windows):
#                print('unique', clones.get_clone(proc))
#
#            print("begin")
#            clones.run_until_complete()
#            print(clones.ln_prob().values())

#            for x, _ in enumerate(windows):
#                jobs.append(mp.Process(target=clones.get_clone(x).run_until_complete()))

            jobs = list()
            criterium = list()
            for _, winx in enumerate(windows):
                criterium.append(fh.criteria_flathist(macro_min=winx[0], macro_max=winx[1]))
            for x, crit in enumerate(criterium):
                jobs.append(mp.Process(target=fh.monte_carlo, args=(x, crit, )))

            for p in jobs:
                p.start()

            for p in jobs:
                p.join()

#            print(clones.ln_prob().values())

if __name__ == "__main__":
    unittest.main()
