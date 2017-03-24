/**
 * \mainpage
 *
 * Developed by Harold Wickes Hatch, 12/13/2013, hhatch.com, harold@hhatch.com
 *
 * An example simulation of a Lennard Jones fluid at T*=1.5
 * Using Transition Matrix Monte Carlo
 *
 * See the NIST Standard Reference Simulation Website
 * http://www.nist.gov/mml/csd/informatics_research/srsw.cfm
 * 
 */

#include "pair_lj.h"
#include "mc_wltmmc.h"

int main(int argc, char** argv) {

  int nMolMax = 370, nMolMin = 0, nfreq = 1e4, ncfreq = 1e6;
  long long int npr = 8000000000;
  double rCut = 3., beta = 1./1.5, lnz = -1.568214, boxl = pow(512, 1./3.);
  std::stringstream addMolType;
  addMolType << "data.lj";

  // parse command-line arguments using getopt
  { int index, c; opterr = 0;
    while ((c = getopt(argc, argv, "b:r:z:l:x:n:f:c:p:m:")) != -1) {
      switch (c) {
        case 'b': beta = atof(optarg); break;
        case 'r': rCut = atof(optarg); break;
        case 'z': lnz = atof(optarg); break;
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
    cout << "# -b beta " << beta << " -r rCut " << rCut << " -z lnz " << lnz << " -l boxl " << boxl << " -x nMolMax " << nMolMax << " -n nMolMin " << nMolMin << " -f nfreq " << nfreq << " -c ncfreq " << ncfreq << " -p npr " << npr << " -m addMolType " << addMolType.str() << endl;
    for (index = optind; index < argc; index++) printf("Non-option argument %s\n", argv[index]);
  }

  // initialize domain, pair interactions and acceptance criteria
  myRanInitByDate();
  Space s(3, 0);
  for (int dim=0; dim < s.dimen(); ++dim) s.lset(boxl,dim);
  stringstream admss;
  admss << s.install_dir() << "/forcefield/" << addMolType.str();
  s.addMolInit(admss.str().c_str());
  PairLJ p(&s, rCut);
  p.initEnergy();
  CriteriaWLTMMC c(beta, exp(lnz),"nmol",nMolMin-0.5,nMolMax+0.5,nMolMax-nMolMin+1);
  
  // initialize MC simulation object
  WLTMMC mc(&s, &p, &c);
  mc.weight=3./4.;
  mc.transformTrial("translate");
  mc.weight=1./8.;
  mc.deleteTrial();
  mc.addTrial(admss.str().c_str());
 
  // output log, lnpi and movie
  mc.initLog("log", nfreq);
  mc.initColMat("colMat", ncfreq); 
  mc.setNFreqCheckE(ncfreq, 2e-4);
  mc.setNFreqTune(ncfreq);
  mc.initMovie("movie", nfreq); 
  //mc.initXTC("movie", nfreq); 
  mc.initRestart("tmp/rst", ncfreq);

  //production tmmc simulation
  c.collectInit();
  c.tmmcInit();
  // nOverlap==-1 makes all OMP windows span entire order parameter range
  // nOverlap>=0 makes OMP windows divide-up range of order parameter 
  // nExp determines the relative size of the windows as a function of processor number
  const int nOverlap = 0;
  const double nExp = 1.5;
  //mc.initWindows(0, -1);
  mc.initWindows(nExp, nOverlap);
  //mc.runNumSweeps(100, -1);
  //mc.runNumTrials(npr);
  mc.runNumSweeps(10, -1); 
  
  cout << "# MC " << s.id() << " elapsed time: " << double(clock()) / double(CLOCKS_PER_SEC) << " seconds" << endl;
}


