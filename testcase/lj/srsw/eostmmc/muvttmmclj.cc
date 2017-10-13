/**
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#include "pair_lj.h"
#include "mc_wltmmc.h"
#include "trial_add.h"
#include "trial_delete.h"
#include "trial_transform.h"
#include "trial_confswap_omp.h"

int main(int argc, char** argv) {

  int openMP = 0;
  int nMolMax = 370, nMolMin = 0, nfreq = 1e4, ncfreq = 1e6;
  int64_t npr = 8000000000;
  double rCut = 3., temp = 1.5, lnz = -1.568214, boxl = 8.;
  std::stringstream molType;
  molType << "data.lj";

  // parse command-line arguments using getopt
  { int index, c; opterr = 0;
    while ((c = getopt(argc, argv, "t:r:z:l:x:n:f:c:p:m:o:")) != -1) {
      switch (c) {
        case 't': temp = atof(optarg); break;
        case 'r': rCut = atof(optarg); break;
        case 'z': lnz = atof(optarg); break;
        case 'l': boxl = atof(optarg); break;
        case 'x': nMolMax = atoi(optarg); break;
        case 'n': nMolMin = atoi(optarg); break;
        case 'f': nfreq = atoi(optarg); break;
        case 'c': ncfreq = atoi(optarg); break;
        case 'o': openMP = atoi(optarg); break;
        case 'p': npr = atoll(optarg); break;
        case 'm': molType.str(""); molType << optarg; break;
	case '?':
          if (isprint(optopt))
            fprintf(stderr, "Unknown option `-%c'.\n", optopt);
          else
            fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
          return 1;
        default: abort();
      }
    }
    cout << "# -t temp " << temp << " -r rCut " << rCut << " -z lnz " << lnz << " -l boxl " << boxl << " -x nMolMax " << nMolMax << " -n nMolMin " << nMolMin << " -f nfreq " << nfreq << " -c ncfreq " << ncfreq << " -p npr " << npr << " -m molType " << molType.str() << " -o openMP " << openMP << endl;
    for (index = optind; index < argc; index++) printf("Non-option argument %s\n", argv[index]);
  }

  // initialize simulation domain
  feasst::ranInitByDate();
  feasst::Space s(3, 0);
  for (int dim=0; dim < s.dimen(); ++dim) s.lset(boxl,dim);
  stringstream addMolType;
  addMolType << s.install_dir() << "/forcefield/" << molType.str().c_str();
  s.addMolInit(addMolType.str().c_str());

  // initialize pair-wise interactions
  feasst::PairLJ p(&s, rCut);
  p.initEnergy();

  // acceptance criteria
  feasst::CriteriaWLTMMC c(1./temp, exp(lnz), "nmol" , nMolMin - 0.5,
                           nMolMax + 0.5, nMolMax - nMolMin + 1);
  c.collectInit();
  c.tmmcInit();

  // initialize MC simulation object
  feasst::WLTMMC mc(&s, &p, &c);
  mc.weight=3./4.;
  feasst::transformTrial(&mc, "translate");
  mc.weight=1./8.;
  feasst::deleteTrial(&mc);
  mc.weight=1./8.;
  feasst::addTrial(&mc, addMolType.str().c_str());

  // if using parallelization, allow configuration swaps between processors
  if (openMP != 0) {
    mc.weight = 1./static_cast<double>(nMolMax);
    mc.confSwapTrial();
  }

  // output log, lnpi and movie
  mc.initLog("log", nfreq);
  mc.initColMat("colMat", ncfreq);
  mc.setNFreqCheckE(ncfreq, 2e-4);
  mc.setNFreqTune(nfreq);
  mc.initMovie("movie", nfreq);
  //mc.initXTC("movie", nfreq);
  mc.initRestart("tmp/rst", ncfreq);

  //production tmmc simulation
  if (openMP == 0) {
    mc.runNumTrials(npr);
  } else {
    mc.initWindows(2.,   // exponent that determines size of windows
                   0);   // extra macrostate overlap between processors
    mc.runNumSweeps(1,   // number of "sweeps"
                   -1);  // maximum number of trials. Infinite if "-1".
  }
}


