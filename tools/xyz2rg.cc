/**
 * \mainpage
 *
 * Developed by Harold Wickes Hatch, 12/13/2013, hhatch.com, harold@hhatch.com
 *
 * This program splits the movie files into separate files corresponding to each bin in the order
 * parameter of the simulation. It also serves as a template for doing some analysis (see
 * pData and dataAccVec). This example trivially computes the average x of first particle.
 *
 * This program reads the log and movie files and assumes that they were generated simultaneously.
 * Make sure that a previous simulation did not start writing a log file, and then the current
 * simulate appended onto this log file, which would cause a mis-match between log and movie files
 *
 * Note that "off by 1" errors are very easy to make, and these analysis codes should be checked
 * that an xyz configuration corresponds with the correct order parameter
 */

#include "pair_tabular.h"
#include "mc_wltmmc.h"

int main(int argc, char** argv) {

  // set input variables
 
  // nMol specifies the desired number of molecules
  //  this is useful in a grand canonical simulation, but default -1 ignores this
  int nMol = -1;

  // skip is used to skip over (skip-1) configurations between each analysis
  int skip = 1;
 
  // these are the default file names
  stringstream ssFileIn("movie"), ssFileOut("analysis"), ssLogFile("log");
 
  // parse command-line arguments using getopt
  { int index, c; opterr = 0;
    while ((c = getopt(argc, argv, "i:o:x:l:")) != -1) {
      switch (c) {
        case 'x': nMol = atoi(optarg); break;
        case 'i': ssFileIn.str(""); ssFileIn << optarg; break;
        case 'o': ssFileOut.str(""); ssFileOut << optarg; break;
        case 'l': ssLogFile.str(""); ssLogFile << optarg; break;
        case '?':
          if (isprint(optopt))
            fprintf(stderr, "Unknown option `-%c'.\n", optopt);
          else
            fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
          return 1;
        default: abort();
      }
    }
    cout << "# -i " << ssFileIn.str() << " -o " << ssFileOut.str() << " -x " << nMol << " -l " << ssLogFile.str() << endl;
    for (index = optind; index < argc; index++) printf("Non-option argument %s\n", argv[index]);
  }

  // initialize simulation
  WLTMMC mc("tmp/rst");
  Space* space = mc.space();
  Pair* pair = mc.pair();
  CriteriaWLTMMC* crit = mc.c();

  // Accumulate averages over all processors
  AccumulatorVec dataAccVec;

  // loop through processors
  const int nProc = mc.nWindows();
  for (int p = 0; p < nProc; ++p) {
   
    // initialize files
    stringstream ss;
    if (nProc == 1) {
      ss << ssFileIn.str() << ".xyz";
    } else {
      ss << ssFileIn.str() << "p" << p << "pr.xyz";
    }
    std::ifstream inFile(ss.str().c_str());
    ss.str("");
    if (nProc == 1) {
      ss << ssLogFile.str();
    } else {
      ss << ssLogFile.str() << "p" << p << "pr";
    }
    std::ifstream logFile(ss.str().c_str());
    ss.str("");
    if (nProc == 1) {
      ss << ssFileOut.str();
    } else {
      ss << ssFileOut.str() << "p" << p << "pr";
    }
    std::ofstream outFile(ss.str().c_str());
 
    //AccumulatorVec pData;
    //vector<bool> firstPrint(crit->nBin(), true);
   
    // prep log file for reading
    if (mc.nFreqMovie() % mc.nFreqLog() != 0) {
      cout << "err: " << mc.nFreqLog() << " " << mc.nFreqMovie()<< endl;
      exit(0);
    }
    const int logSkip = myRound(mc.nFreqMovie() / mc.nFreqLog());
    string line;
    getline(logFile, line);
  
    // read xyz
    int iter = 0;
    while (!inFile.eof()) {
      Space* spacetmp = space->clone();
      Pair* pairtmp = pair->clone(spacetmp);
      pairtmp->readxyz(inFile);
     
      // skip log file lines
      for (int i = 0; i < logSkip ; ++i) getline(logFile, line);
     
      // check for EOF and number of molecules
      if ( (spacetmp->natom() != 0) && ( (nMol == -1) || (spacetmp->nMol() == nMol) ) ) {
        if ( (!inFile.eof()) && (iter%skip==0) ) {
         
          // parse log file lines
          //std::istringstream iss(line);
          //vector<double> data(5);
          //iss >> data[0] >> data[1] >> data[2] >> data[3] >> data[4];
          //const double beta = data[3];
          //if ( (beta < crit->mMax()) && (beta > crit->mMin()) ) { //catch out of bounds due to EOF in log
          //  const int bin = crit->bin(beta);
          const int bin = spacetmp->nMol();
          if ( (bin < crit->mMax()) && (bin > crit->mMin()) ) { //catch out of bounds due to EOF in log

            // compute the radius of gyration for a given bin
            for (int iMol = 0; iMol < bin; ++iMol) {
              double rg = 0;
              const int firstAtom = spacetmp->mol2part()[iMol];
              int sites = 0;
              for (int iAtom = firstAtom+1; iAtom < spacetmp->mol2part()[iMol+1]; ++iAtom) {
                ++sites;
                for (int dim = 0; dim < spacetmp->dimen(); ++dim) {
                  rg += pow(spacetmp->x(iAtom, dim)-spacetmp->x(firstAtom, dim), 2.);
                }
              }
              cout << "rg " << rg/double(sites) << " bin " << bin << endl;
              dataAccVec.accumulate(bin, rg/double(sites));
            }
          }
        }
      }
      spacetmp->cellOff();  // avoid error check
    }
   
//    // output processor specific averages
//    outFile << "# " << crit->mType() << endl;
//    for (int i = 0; i < dataAccVec.size(); ++i) {
//      outFile << crit->bin2m(i) << " " << dataAccVec.vec()[i].average() << endl;
//    }
//
//    // accumulate over all processors
//    for (int i = 0; i < pData.size(); ++i) {
//      dataAccVec.accumulate(i, pData.vec(i).average());
//    }
  }

  // output average over all p
  stringstream ss;
  ss << ssFileOut.str();
  std::ofstream outFile(ss.str().c_str());
  outFile << "# " << crit->mType() << "stdev" << endl;
  for (int i = 0; i < dataAccVec.size(); ++i) {
    outFile << crit->bin2m(i) << " " << dataAccVec.vec(i).average() << " "<< dataAccVec.vec(i).stdev() << endl;
  }
}
