/**
 * \mainpage
 *
 * Developed by Harold Wickes Hatch, 12/13/2013, hhatch.com, harold@hhatch.com
 *
 *
 */

#include "pair_lj.h"
#include "mc_wltmmc.h"
#include "trial_add.h"
#include "trial_delete.h"
#include "trial_transform.h"

using namespace feasst;

int main(int argc, char** argv) {

  int nMolMax = 370, nMolMin = 0, nfreq = 1e4, ncfreq = 1e6;
  int64_t npr = 8000000000;
  double rCut = 3., beta = 1./0.7, activ = exp(-5.943376), boxl = pow(512, 1./3.);
  std::stringstream addMolType;
  addMolType << "data.lj";

  // parse command-line arguments using getopt
  { int index, c; opterr = 0;
    while ((c = getopt(argc, argv, "b:r:a:l:x:n:f:c:p:m:")) != -1) {
      switch (c) {
        case 'b': beta = atof(optarg); break;
        case 'r': rCut = atof(optarg); break;
        case 'a': activ = atof(optarg); break;
        case 'l': boxl = atof(optarg); break;
        case 'x': nMolMax = atoi(optarg); break;
        case 'n': nMolMin = atoi(optarg); break;
        case 'f': nfreq = atoi(optarg); break;
        case 'c': ncfreq = atoi(optarg); break;
        case 'p': npr = atoll(optarg); break;
        case 'm': addMolType.str(""); addMolType << optarg; break;
	case '?':
          if (isprint(optopt))
            fprintf(stderr, "Unknown option `-%c'.\n", optopt);
          else
            fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
          return 1;
        default: abort();
      }
    }
    cout << "# -b beta " << beta << " -r rCut " << rCut << " -a activ " << activ << " -l boxl " << boxl << " -x nMolMax " << nMolMax << " -n nMolMin " << nMolMin << " -f nfreq " << nfreq << " -c ncfreq " << ncfreq << " -p npr " << npr << " -m addMolType " << addMolType.str() << endl;
    for (index = optind; index < argc; index++) printf("Non-option argument %s\n", argv[index]);
  }

  // initialize domain, pair interactions and acceptance criteria
  ranInitByDate();
  Space s(3, 0);
  for (int dim=0; dim < s.dimen(); ++dim) s.lset(boxl,dim);
  stringstream addMolSS;
  addMolSS << s.install_dir() << "/forcefield/"<< addMolType.str();
  s.addMolInit(addMolSS.str().c_str());
  PairLJ p(&s, rCut);
  p.initEnergy();
  CriteriaWLTMMC c(beta, activ,"nmol",nMolMin-0.5,nMolMax+0.5,nMolMax-nMolMin+1);

  // initialize MC simulation object
  WLTMMC mc(&s, &p, &c);
  mc.weight=3./4.;
  transformTrial(&mc, "translate");
  mc.weight=1./8.;

  //deleteTrial(&mc);
  //mc.initTrial(new TrialDelete(addMolSS.str().c_str()));
  //mc.initTrial(new TrialAdd(addMolSS.str().c_str()));
  deleteTrial(&mc);
  addTrial(&mc, addMolSS.str().c_str());

  // output log, lnpi and movie
  mc.initLog("log", nfreq);
  mc.initColMat("colMat", ncfreq);
  mc.setNFreqCheckE(ncfreq, 2e-4);
  mc.setNFreqTune(nfreq);
  mc.initMovie("movie", nfreq);
  //mc.initXTC("movie", nfreq);
  mc.initRestart("tmp/rst", ncfreq);

  //production tmmc simulation
  c.collectInit();
  c.tmmcInit();
  //mc.runNumTrials(npr);
  mc.initWindows(1);
  //mc.nMolSeekInRange();
  mc.runNumSweeps(10, -1);

  cout << "# MC " << s.id() << " elapsed time: " << double(clock()) / double(CLOCKS_PER_SEC) << " seconds" << endl;
}


