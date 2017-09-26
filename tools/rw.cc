/**
 * Reweight a collection matrix file to a given lnz, if given.
 * If a lnz is not given, then attempt to find saturation.
 *   For saturation, if a volume is not given, attempt to find the volume from a
 *   restart file in order to compute the density and pressure at saturation.
 */

#include "functions.h"
#include "pair_ideal.h"
#include "mc_wltmmc.h"

int main(int argc, char** argv) {

  // set input variables
  double volume = -1;
  const double LNZDEFAULT = 11234533;
  double lnz = LNZDEFAULT;
  int phaseBoundary = -1;
  std::ostringstream rstFile("tmp/rstspace"), ssFileIn("colMat.txt"), ssFileOut("colMatrw.txt");

  // parse command-line arguments using getopt
  { int index, c; opterr = 0;
    while ((c = getopt(argc, argv, "v:z:p:r:i:o:s:")) != -1) {
      switch (c) {
        case 'v': volume = atof(optarg); break;
        case 'z': lnz = atof(optarg); break;
        case 'p': phaseBoundary = atoi(optarg); break;
        case 'r': rstFile.str(""); rstFile << optarg; break;
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
    cout << "#"
      << " -i " << ssFileIn.str()
      << " -o " << ssFileOut.str()
      << " -z " << lnz
      << " -p " << phaseBoundary
      << " -v " << volume
      << " -r " << rstFile.str()
      << endl;
    for (index = optind; index < argc; index++) printf("Non-option argument %s\n", argv[index]);
  }

  // If the lnz was set to some value, reweight to the given lnz.
  if (lnz != LNZDEFAULT) {
    feasst::CriteriaWLTMMC criteria(ssFileIn.str().c_str());
    criteria.readCollectMat(ssFileIn.str().c_str());
    criteria.lnPIrw(exp(lnz));
    criteria.printRWinit();
    criteria.printCollectMat(ssFileOut.str().c_str());

  // Otherwise, attempt to find saturation and reweight the lnzsat.
  } else {
    // if volume is not given or unphysical, attempt to find volume from restart
    shared_ptr<feasst::Space> space;
    if (volume < 0) {
      ASSERT(feasst::fileExists(rstFile.str().c_str()), "if volume is not "
        << "specified, then a restart file must be provided");
      space = make_shared<feasst::Space>(rstFile.str().c_str());
    } else {
      space = make_shared<feasst::Space>();
      space->lset(1);
      space->scaleDomain(volume);
    }

    feasst::PairIdeal pair(space.get(), 0.);
    feasst::CriteriaWLTMMC criteria(ssFileIn.str().c_str());
    feasst::WLTMMC wltmmc(space.get(), &pair, &criteria);
    criteria.readCollectMat(ssFileIn.str().c_str());
    if (phaseBoundary != -1) {
      criteria.setPhaseBoundary(phaseBoundary);
    }
    wltmmc.printSat();
  }






//  feasst::WLTMMC mc("tmp/rst");
//  mc.c()->readCollectMat(ssFileIn.str().c_str());
//  mc.initLog("logtest", 0);
//  mc.c()->lnPIrw(exp(lnz));
//  mc.c()->printRWinit();
//  mc.c()->printCollectMat(ssFileOut.str().c_str());

//  mc.printSat();
}
