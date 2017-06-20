/**
 * \mainpage
 *
 * Developed by Harold Wickes Hatch, 12/13/2013, hhatch.com, harold@hhatch.com
 *
 *
 */

#include "space.h"

int main(int argc, char** argv) {

  // set input variables
  int nMol = 500, preMicellarAgg = 5;
  stringstream ssFileIn("test/prod.xtc"), ssFileOut("cluster.txt"), molType("/home/hwh/hwhcms/forcefield/data.cg3_60_43_1alt");
 
  // parse command-line arguments using getopt
  { int index, c; opterr = 0;
    while ((c = getopt(argc, argv, "p:n:i:o:m:")) != -1) {
      switch (c) {
        case 'i': ssFileIn.str(""); ssFileIn << optarg; break;
        case 'o': ssFileOut.str(""); ssFileOut << optarg; break;
        case 'm': molType.str(""); molType << optarg; break;
        case 'n': nMol = atoi(optarg); break;
        case 'p': preMicellarAgg = atoi(optarg); break;
        case '?':
          if (isprint(optopt))
            fprintf(stderr, "Unknown option `-%c'.\n", optopt);
          else
            fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
          return 1;
        default: abort();
      }
    }
    cout << "# -i " << ssFileIn.str() << " -o " << ssFileOut.str() << " -n " << nMol << " -m " << molType.str() << " -p " << preMicellarAgg << endl;
    for (index = optind; index < argc; index++) printf("Non-option argument %s\n", argv[index]);
  }

  // initialize xtc input file
  XDRFILE* trjFileXDR;
  trjFileXDR = xdrfile_open(ssFileIn.str().c_str(), "r");

  // initialize output
  std::ofstream outFile(ssFileOut.str().c_str());
  outFile << "# natom nClusters avClusterSize asphericity acylindricity relShpAniso Rg" << endl;
 
  // initialize number of molecules, read first frame xtc
  Space s("tmp/rstspace");
  for (int i = 0; i < nMol; ++i) s.addMol(0);
  int endXTC = s.readXTC(ssFileIn.str().c_str(), trjFileXDR);
 
  // initialize cell list and cluster definition
  s.initCellAtomCut(1);
  s.updateCells(4./3.);
  s.addTypeForCluster(1);
  s.preMicellarAgg(preMicellarAgg);

  // read coordinates from xtc until eof
  while (endXTC == 0) {
     
    // compute clusters
    s.wrapMol();
    s.updateCellofallMol();
    s.updateClusters(4./3.);
    s.xClusterGen();
    s.xClusterShape();
    outFile << s.nMol() << " " << s.nClusters() << " " << s.clusterAvSize() << " " << s.clusterAsphericityAv() << " " << s.clusterAcylindricityAv() << " " << s.clusterRelShapeAnisoAv() << " " << s.clusterRgAv() << endl;

    // attempt to read next trajectory frame
    endXTC = s.readXTC(ssFileIn.str().c_str(), trjFileXDR);
  }

  // output final cluster statistics
  outFile << "# avClusterNum avClusterSize freemon" << endl;
  outFile << s.clusterNumAccVec().vec(s.nMol()).average() << " " << s.clusterSizeAccVec().vec(s.nMol()).average() << " " << s.freeMon().average() << endl;
 
  stringstream ssFileOut2;
  ssFileOut2 << ssFileOut.str() << "2";
  s.printClusterStat(ssFileOut2.str().c_str());
 
  // clean up
  s.cellOff();  // avoid error check
  xdrfile_close(trjFileXDR);
}


