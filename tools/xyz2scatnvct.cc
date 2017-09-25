/**
 * \mainpage
 *
 * Developed by Harold Wickes Hatch, 12/13/2013, hhatch.com, harold@hhatch.com
 *
 *
 */

#include "pair_tabular.h"
#include "mc_wltmmc.h"
#include "analyze_scatter.h"

int main(int argc, char** argv) {

  // set input variables
  const double np = 12;
  double boxl = 9, rCut = 4./3.;
  int preMicellarAgg = 5, nMolMax = -1, skip = 100;
  stringstream ssFileIn("movie"), ssFileOut("cluster.txt"), ssLogFile("log");
 
  // parse command-line arguments using getopt
  { int index, c; opterr = 0;
    while ((c = getopt(argc, argv, "r:i:o:b:p:x:")) != -1) {
      switch (c) {
        case 'b': boxl = atof(optarg); break;
        case 'r': rCut = atof(optarg); break;
        case 'p': preMicellarAgg = atoi(optarg); break;
        case 'x': nMolMax = atoi(optarg); break;
        case 'i': ssFileIn.str(""); ssFileIn << optarg; break;
        case 'o': ssFileOut.str(""); ssFileOut << optarg; break;
        case '?':
          if (isprint(optopt))
            fprintf(stderr, "Unknown option `-%c'.\n", optopt);
          else
            fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
          return 1;
        default: abort();
      }
    }
    cout << "# -i " << ssFileIn.str() << " -o " << ssFileOut.str() << " -b " << boxl << " -r " << rCut << " -p " << preMicellarAgg << " -x " << nMolMax << endl;
    for (index = optind; index < argc; index++) printf("Non-option argument %s\n", argv[index]);
  }

  // initialize simulation
  MC mc("tmp/rst");
  Space* space = mc.space();
  Pair* pair = mc.pair();
  space->initEuler(1);

  // accumulate averages over all p
  AccumulatorVec clusterphi, coordphi;
 
  // initialize files
  stringstream ss;
  ss << "movie.xyz";
  std::ifstream inFile(ss.str().c_str());
  ss.str("");
  ss << "scatpost";

  AnalyzeScatter scat(space, pair);
  superball sball(0.5, 0.5, 0.5, 1., 1.);
  sball.genBeadRep(0.1, 5);
  scat.initAniso(&sball);
  scat.initFileName(ss.str().c_str());
  scat.initSANS(0.0001);

  // read xyz
  int iter = 0;
  while (!inFile.eof()) {
    pair->readxyz(inFile);
    if ( (space->natom() != 0) && ( (nMolMax == -1) || (space->nMol() == nMolMax) ) ) {
      if ( (!inFile.eof()) && (iter%skip==0) ) {
       
        // compute analysis
        scat.update();
      }
    }
    ++iter;
    //cout << "iter " << iter << endl;
  }
  scat.print();
  space->cellOff();  // avoid error check
}
