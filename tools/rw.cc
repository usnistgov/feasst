/**
 * \mainpage
 *
 * Developed by Harold Wickes Hatch, 12/13/2013, hhatch.com, harold@hhatch.com
 *
 *
 */

#include "mc_wltmmc.h"

int main(int argc, char** argv) {

  // set input variables
  double temp = 0, activ = 0;
  int nMolMax = -1, nMolMin = 0, nSmooth = 1;
  temp = 525; activ = exp(-8.08564); nMolMax = 265; nMolMin = 0;
  std::ostringstream molType("a"), ssFileIn("colMat.txt"), ssFileOut("colMatrw.txt");
 
  // parse command-line arguments using getopt
  { int index, c; opterr = 0;
    while ((c = getopt(argc, argv, "t:a:x:n:m:i:o:s:")) != -1) {
      switch (c) {
        case 't': temp = atof(optarg); break;
        case 'a': activ = exp(atof(optarg)); break;
        case 'x': nMolMax = atoi(optarg); break;
        case 'n': nMolMin = atoi(optarg); break;
        case 's': nSmooth = atoi(optarg); break;
        case 'm': molType.str(""); molType << optarg; break;
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
    cout << "# -t " << temp << " -a " << activ << " -x " << nMolMax << " -n " << nMolMin << " -i " << ssFileIn.str() << " -o " << ssFileOut.str() << " -m " << molType.str() << " -s " << nSmooth << endl;
    for (index = optind; index < argc; index++) printf("Non-option argument %s\n", argv[index]);
  }
 
  WLTMMC mc("tmp/rst");
  mc.c()->readCollectMat(ssFileIn.str().c_str());
  mc.initLog("logtest", 0);
  mc.c()->lnPIrw(activ);
  mc.c()->printRWinit();
  mc.c()->printCollectMat(ssFileOut.str().c_str());

//  mc.printSat();
}
