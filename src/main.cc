/**
 * \mainpage
 *
 * Developed by Harold Wickes Hatch, 12/13/2013, hhatch.com, harold@hhatch.com
 *
 * 
 */

#include "assert.h"
#include "pair_lj_multi.h"
#include "mc_wltmmc.h"
#include "ui_abbreviated.h"
#include "trial_swap.h"
#include "analyze.h"

class AnalyzeComp : public Analyze {
 public:
  AnalyzeComp(Space *space, Pair *pair) : Analyze(space, pair) {}
  ~AnalyzeComp() {}
  Accumulator nA;
  void update() {
  //void update(const int iMacro) {
    nA.accumulate(space_->nMolType()[0]);
  }
};

int main(int argc, char** argv) {

  int nMolMax = 38, nMolMin = 0, nfreq = 1e4, ncfreq = 1e5;
  int64_t npr = 1e6;
  double rCut = 3., beta = 1./1.5, activ = exp(-4.568214), boxl = pow(512, 1./3.);
  std::stringstream addMolTypeA, addMolTypeB;
  Space stmp(3, 0);
  addMolTypeA << stmp.install_dir() << "/forcefield/data.lj";
  addMolTypeB << stmp.install_dir() << "/forcefield/data.ljb";

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
	case '?':
          if (isprint(optopt))
            fprintf(stderr, "Unknown option `-%c'.\n", optopt);
          else
            fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
          return 1;
        default: abort();
      }
    }
    cout << "# -b beta " << beta << " -r rCut " << rCut << " -a activ " << activ << " -l boxl " << boxl << " -x nMolMax " << nMolMax << " -n nMolMin " << nMolMin << " -f nfreq " << nfreq << " -c ncfreq " << ncfreq << " -p npr " << npr << endl;
    for (index = optind; index < argc; index++) printf("Non-option argument %s\n", argv[index]);
  }

  // initialize domain, pair interactions and acceptance criteria
  myRanInitByDate();
  myRanInitForRepro(14902901300);
//  myRanInitForRepro(1490289674);
  Space s(3, 0);
  for (int dim=0; dim < s.dimen(); ++dim) s.lset(boxl,dim);
  PairLJMulti p(&s, rCut);
  s.addMolInit(addMolTypeA.str().c_str());
  p.initLMPData(addMolTypeA.str().c_str());
  s.addMolInit(addMolTypeB.str().c_str());
  p.initLMPData(addMolTypeB.str().c_str());
  p.setLambdaij(0, 0, 1);
  p.setLambdaij(0, 1, 0.5);
  p.setLambdaij(1, 1, 1);
  p.lrcFlag = 0;
  p.cutShift(1);
  p.initEnergy();
  
  // initialize acceptance criteria
  CriteriaMetropolis c(beta, activ);
  //c.activVec.push_back(exp(-3.));   //!< activity of B 
  c.addActivity(exp(-3));

  // initialize MC simulation object
  MC mc(&s, &p, &c);
  mc.weight=3./4.;
  transformTrial(&mc, "translate");
  
  // generate initial configuration
  assert(nMolMax % 2 == 0);
  mc.nMolSeek(nMolMax/2., addMolTypeA.str().c_str(), 1e8);
  mc.nMolSeek(nMolMax, addMolTypeB.str().c_str(), 1e8);

  // add more trials
  mc.weight=1./4.;
  shared_ptr<TrialSwap> tswp = make_shared<TrialSwap>(
    addMolTypeA.str().c_str(), addMolTypeB.str().c_str());
  mc.initTrial(tswp);

  // output log, lnpi and movie
  mc.initLog("log", nfreq);
  mc.setNFreqCheckE(ncfreq, 2e-4);
  mc.setNFreqTune(nfreq);
  mc.initMovie("movie", nfreq); 
  //mc.initXTC("movie", nfreq); 
  mc.initRestart("tmp/rst", ncfreq);

  // analysis
  shared_ptr<AnalyzeComp> a = make_shared<AnalyzeComp>(&s, &p);
  mc.initAnalyze(a.get());

  //production tmmc simulation
  mc.runNumTrials(npr);
  
  //cout << "av nA " << tswp->nA().average() << " +/- " << tswp->nA().blockStdev() << endl;
  //cout << "av nB " << tswp->nB().average() << " +/- " << tswp->nB().blockStdev() << endl;
  //cout << "xA " << tswp->nA().average()/s.nMol() << " +/- " << tswp->nA().blockStdev()/s.nMol() << endl;
  cout << "an " << a->nA.average()/double(s.nMol()) << " +/- " << a->nA.blockStdev()/double(s.nMol()) << endl;
//  // alternatively, accumulate the composition
//  Accumulator xA;
//  const int trialPerAccum = 10;
//  for (int step = 0; step < npr/trialPerAccum; ++step) {
//    mc.runNumTrials(trialPerAccum);
//    xA.accumulate(static_cast<double>(s.nMolType()[0])/
//                  static_cast<double>(s.nMolType()[1]));
//  }

  cout << "# MC " << s.id() << " elapsed time: " << double(clock()) / double(CLOCKS_PER_SEC) << " seconds" << endl;
}


