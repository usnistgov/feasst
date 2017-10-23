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

// Compare computed macrostate with SRSW
void compareEnergyAndMacro(feasst::CriteriaWLTMMC criteria,
  int iMacro, //!< macrostate (e.g., number of particles)
  double peAv, double peStd, double lnPIav, double lnPIstd) {
  ASSERT(criteria.mMax() > iMacro, "comparing macrostate(" << iMacro
    << ") that doesn't exist in criteria, with max of " << criteria.mMax());
  const double diff = criteria.pe(iMacro).average() - peAv;
  // 99% confidence interval
  const double tol = 2.576*(criteria.pe(iMacro).blockStdev() + peStd);
  ASSERT(fabs(diff) < tol,
    "N=" << iMacro << " energy is " << criteria.pe(iMacro).average() << " +/- "
    << criteria.pe(iMacro).blockStdev() << " but SRSW is " << peAv << " +/- "
    << peStd << "." << endl << "The difference is: " << diff << endl
    << "The tolerance (99%) is: " << tol);
}

int main(int argc, char** argv) {  // LJ, SRSW_EOSTMMC

  int openMP = 0;
  int nMolMax = 5, nMolMin = 0, nfreq = 1e4, ncfreq = 1e6;
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
    cout << "# -t temp " << temp << " -r rCut " << rCut << " -z lnz " << lnz << " -l boxl " << boxl << " -x nMolMax " << nMolMax << " -n nMolMin " << nMolMin << " -f nfreq " << nfreq << " -c ncfreq " << ncfreq << " -m molType " << molType.str() << " -o openMP " << openMP << endl;
    for (index = optind; index < argc; index++) printf("Non-option argument %s\n", argv[index]);
  }  // GETOPT

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

  // output log, lnpi and movie
  mc.initLog("log", nfreq);
  mc.initColMat("colMat", ncfreq);
  mc.setNFreqCheckE(ncfreq, 2e-4);
  mc.setNFreqTune(nfreq);
  mc.initMovie("movie", nfreq);
  //mc.initXTC("movie", nfreq);
  mc.initRestart("tmp/rst", ncfreq);

  //production tmmc simulation
  if (openMP != 0) {
    mc.initWindows(2.,   // exponent that determines size of windows
                   0);   // extra macrostate overlap between processors
  }
  mc.runNumSweeps(20,   // number of "sweeps"
                 -1);  // maximum number of trials. Infinite if "-1".

  // Test results against the SRSW values
  if (openMP == 0) {
    ASSERT(fabs(c.pe(0).average()) < 1e-13,
      "N=0 should have zero energy, but pe=" << c.pe(0).average());

    compareEnergyAndMacro(c, 1,
      -0.0006057402333333332,
      6.709197666659334e-10,
      -270.0061768+274.6763737666667,
      0.037092307087640365);

    compareEnergyAndMacro(c, 2,
      -0.030574223333333334,
      9.649146611661053e-06,
      -266.0191155333334+274.6763737666667,
      0.03696447428346385);

    compareEnergyAndMacro(c, 3,
      -0.089928316,
      0.0001387472078025413,
      -262.4277240666667+274.6763737666667,
      0.037746391500313385);

    compareEnergyAndMacro(c, 4,
      -0.1784570533333333,
      3.3152449884326804e-05,
      -259.11444086666665+274.6763737666667,
      0.03809721387875822);

    compareEnergyAndMacro(c, 5,
      -0.29619201333333334,
      1.3487910636322294e-05,
      -256.0144809+274.6763737666667,
      0.03845757460933292);
  }
}


