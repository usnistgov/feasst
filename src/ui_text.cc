#include "./ui_abbreviated.h"
#include "./pair_lj_multi.h"

#ifdef FEASST_NAMESPACE_
using namespace feasst;
#endif  // FEASST_NAMESPACE_

shared_ptr<Space> space;
shared_ptr<Criteria> criteria;
shared_ptr<Pair> pair;
shared_ptr<MC> mc;
shared_ptr<WLTMMC> wltmmc;

void interpret(vector<string> cmd) {
  if (cmd[0] == "space") {
    const int dimen = stoi(cmd[1]);
    space = make_shared<Space>(dimen, 0);

  } else if (cmd[0] == "boxl") {
    double boxl = stod(cmd[1]);
    for (int dim = 0; dim < space->dimen(); ++dim) {
      space->lset(boxl, dim);
    }

  } else if (cmd[0] == "addMolInit") {
    space->addMolInit(cmd[1].c_str());

  } else if (cmd[0] == "pair") {
    if (cmd[1] == "lj") {
      const double rCut = stod(cmd[2]);
      pair = make_shared<PairLJ>(space.get(), rCut);
    }

  } else if (cmd[0] == "initEnergy") {
    pair->initEnergy();

  } else if (cmd[0] == "criteria") {
    const double beta = stod(cmd[2]);
    const double lnz = stod(cmd[3]);
    if (cmd[1] == "metropolis") {
      criteria = make_shared<CriteriaMetropolis>(beta, lnz);
    }

  } else if (cmd[0] == "mc") {
    mc = make_shared<MC>(space.get(), pair.get(), criteria.get());

  } else if (cmd[0] == "trial") {
    if (cmd[1] == "translate") {
      const double maxMoveParam = stod(cmd[2]);
      transformTrial(mc, cmd[1].c_str(), maxMoveParam);
    }

  } else if (cmd[0] == "nMolSeek") {
    const int nMol = stoi(cmd[1]);
    mc->nMolSeek(nMol, cmd[2].c_str());

  } else if (cmd[0] == "log") {
    mc->initLog(cmd[1].c_str(), stoi(cmd[2]));

  } else if (cmd[0] == "movie") {
    mc->initMovie(cmd[1].c_str(), stoi(cmd[2]));

  } else if (cmd[0] == "run") {
    const int nTrials = stoi(cmd[1]);
    mc->runNumTrials(nTrials);

  }
}

void interpret(std::string line) {
  std::string delimiter = " ";
  size_t pos = 0;
  string token;
  bool splitting = true;
  int splits = 0;
  vector<string> splitString;
  while (splitting && (splits < 1e5)) {
  //while ((pos = line.find(delimiter)) != std::string::npos) {
    if (!((pos = line.find(delimiter)) != std::string::npos)) {
      splitting = false;
    }
    token = line.substr(0, pos);
    splitString.push_back(token);
//    std::cout << token << std::endl;
    line.erase(0, pos + delimiter.length());
    ++splits;
  }
  ASSERT(splits <= 1e5, "error reading line");
  interpret(splitString);
}

int main(int argc, char** argv) {
  std::stringstream fileName;

  // parse command-line arguments using getopt
  { int index, c; opterr = 0;
    while ((c = getopt(argc, argv, "b:r:a:l:x:n:f:c:p:m:")) != -1) {
      switch (c) {
        case 'f': fileName.str(""); fileName << optarg; break;
        case '?':
          if (isprint(optopt))
            fprintf(stderr, "Unknown option `-%c'.\n", optopt);
          else
            fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
          return 1;
        default: abort();
      }
    }
    cout << "# -f fileName " << fileName.str() << endl;
    for (index = optind; index < argc; index++) printf("Non-option argument %s\n", argv[index]);
  }

  std::ifstream input(fileName.str().c_str());
  std::string line;
  int nline = 1;
//  const char commentChar = "#";

  while (std::getline(input, line)) {
    if (!line.empty()) {
      //if (strcmp(*line.at(0), "#") != 0) {
      if (line.at(0) != '#') {
        cout << "# line" << nline << ": " << line << endl;
        interpret(line);
      }
    }
    ++nline;
  }
  return 1;
}
