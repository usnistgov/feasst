/**
 * \file
 *
 * \brief randomly selects monte carlo trials
 *
 */

#include "./mc_wltmmc.h"

/**
 * Constructor
 */
WLTMMC::WLTMMC(Space* space,
  Pair* pair,
  CriteriaWLTMMC* criteria
  )
    : MC(space, pair, criteria),
      c_(criteria) {
  defaultConstruction();
}
WLTMMC::WLTMMC(const char* fileName)
    : MC(fileName) {
  defaultConstruction();

  // cout << "initialize criteria" << endl;
  string strtmp = fstos("rstFileCriteria", fileName);
  ASSERT(!strtmp.empty(), "file name of criteria restart file not provided");
  c_ = new CriteriaWLTMMC(strtmp.c_str());
  delete criteria_;
  criteria_ = c_;
  criteriaOwned_ = true;

  // inialize flags associated with trials
  strtmp = fstos("densThresConfigBias", fileName);
  if (!strtmp.empty()) {
    densThresConfigBias_ = atoi(strtmp.c_str());
  }

  // cout << "reinitialize criteria in trials" << endl;
  for (int i = 0; i < nTrials(); ++i) {
    trialVec_[i]->replaceCriteria(c_);
  }

  colMatFileName_ = fstos("colMatFileName", fileName);
  nFreqColMat_ = fstoi("nFreqColMat", fileName);

  strtmp = fstos("nWindow", fileName);
  if (!strtmp.empty()) {
    nWindow_ = stoi(strtmp);
    nExp_ = fstod("nExp", fileName);
    nOverlap_ = fstoi("nOverlap", fileName);
    betaInc_ = fstoi("betaInc", fileName);
    lnzInc_ = fstoi("lnzInc", fileName);
  }

  // read radial distribution function restart files
  strtmp = fstos("nFreqGR", fileName);
  if (!strtmp.empty()) {
    nFreqGR_ = stoi(strtmp);
    GRFileName_ = fstos("GRFileName", fileName);
    const int nfiles = fstoi("numGRFiles", fileName);
    gr_.resize(nfiles);
    for (int i = 0; i < nfiles; ++i) {
      for (int n = 0; n < c()->nBin(); ++n) {
        stringstream ss;
        ss << fileName << "gri" << i << "n" << n;
        if (myFileExists(ss.str().c_str())) {
          gr_[i].push_back(make_shared<Histogram>(ss.str().c_str()));
        }
      }
      grt_.push_back(make_shared<Histogram>(gr_[i].front()->binWidth(),
                                            gr_[i].front()->iType(),
                                            gr_[i].front()->jType()));
    }
  }

  strtmp = fstos("wlFlatProd", fileName);
  if (!strtmp.empty()) {
    wlFlatProd_ = stoi(strtmp);
  }
  strtmp = fstos("wlFlatTerm", fileName);
  if (!strtmp.empty()) {
    wlFlatTerm_ = stoi(strtmp);
  }

  if (c_->nSweep() > 1) production_ = true;
}

/**
 * defaults in constructor
 */
void WLTMMC::defaultConstruction() {
  className_.assign("WLTMMC");
  verbose_ = 0;
  nFreqColMat_ = 1e6;
  nFreqGR_ = 0;
  colMatFileName_.assign("colMat");
  initWindows(0);
  betaInc_ = 0;
  lnzInc_ = 0;
  densThresConfigBias_ = 0;
  nMolSeekTarget_ = -1;
  wlFlatProd_ = -1;
  wlFlatTerm_ = -1;
}

WLTMMC::~WLTMMC() {
  // if (pair_ != NULL) pair_->checkEnergy(1e-9, 0);
  checkTrialCriteria();
  if (criteriaOwned_) {
    delete c_;
  }
}

/**
 * reset object pointers
 */
void WLTMMC::reconstruct() {
  CriteriaWLTMMC* c = c_->clone();
  criteriaOwned_ = true;
  c_ = c;  // swap temporaries for exception safety
  criteria_ = reinterpret_cast<Criteria*>(c_);

  // reconstruct the shared pointers to the radial distribution functions
  for (unsigned int itype = 0; itype < gr_.size(); ++itype) {
    grt_[itype] = grt_[itype]->cloneShrPtr();
    for (unsigned int n = 0; n < gr_[itype].size(); ++n) {
      gr_[itype][n] = gr_[itype][n]->cloneShrPtr();
    }
  }

  MC::reconstruct();
}

/**
 * clone design pattern
 */
WLTMMC* WLTMMC::clone() const {
  WLTMMC* mc = new WLTMMC(*this);
  mc->reconstruct();
  return mc;
}

/**
 * clone design pattern
 */
shared_ptr<MC> WLTMMC::cloneImpl() const {
  shared_ptr<WLTMMC> mc = make_shared<WLTMMC>(*this);
  mc->reconstruct();
  return mc;
}

/**
 * this function is called after every trial attempt
 */
void WLTMMC::afterAttempt() {
  if ( ( (c_->nSweep() > 1) ||
         ((wlFlatProd_ != -1) && (c_->wlFlat() >= wlFlatProd_) ) )
    && (!production_) ) {
    initProduction();
  }
  afterAttemptBase();
  if (nAttempts_ % nFreqColMat_ == 0) {
    c_->lnPIupdate();
    if (!colMatFileName_.empty()) {
//      c_->lnPIpressureIso(space_->vol());
      c_->printCollectMat(colMatFileName_.c_str());

//      // print colmats for each step to animate plot
//      stringstream ss;
//      ss << colMatFileName_ << nAttempts_/nFreqColMat_;
//      c_->printCollectMat(ss.str().c_str());
    }

    // print cluster statistics
    if (space_->clusterType().size() != 0) {
      stringstream ss;
      ss << logFileName_ << "clus";
      space_->printClusterStat(ss.str().c_str());
    }
  }

  // compute and print radial distribution function
  if ( (nFreqGR_ != 0) && (production_) ) {
    if (nAttempts_ % nFreqGR_ == 0) {
      if (!GRFileName_.empty()) {
        const int nBin = c()->bin(c()->mOld());
        for (unsigned int i = 0; i < grt_.size(); ++i) {
          while (static_cast<int>(gr_[i].size()) < nBin+1) {
            gr_[i].push_back(grt_[i]->cloneShrPtr());
          }
          space_->nRadialHist(gr_[i][nBin].get());
          stringstream ss;
          ss << GRFileName_ << "i" << i << "m" << c()->bin2m(nBin);
          space_->printRadial(*gr_[i][nBin], ss.str().c_str());
        }
      }
    }
  }

  // run analyzers
  for (vector<Analyze*>::iterator it = analyzeVec_.begin();
       it != analyzeVec_.end(); ++it) {
    if (nAttempts_ % (*it)->nFreq() == 0) {
      (*it)->update(c_->iMacro());
    }
    if (nAttempts_ % (*it)->nFreqPrint() == 0) {
      (*it)->print(c_);
    }
  }
}

/**
 * determine maximum number of particles for a given temperature and large activity
 */
int WLTMMC::nMolMax(const long long npr,   //!< number of steps in production run
  const double activ,       //!< large activity
  const int nMolExtra       //!< add extra mols
  ) {
  const double mMax = (c_->mMax() - 0.5) * 2 + 0.5;
  c_->initBins(mMax, c_->mMin(), mMax -c_->mMin());
  const int nmx = MC::nMolMax(npr, activ, nMolExtra);
  const double mMax2 = nmx+0.5;
  c_->initBins(mMax2, c_->mMin(), mMax2 -c_->mMin());
  return nmx;
}

/**
 * run for a number of sweeps
 */
void WLTMMC::runNumSweeps(const int nSweeps,  //!< target number of sweeps
  const long long nprMax  //!< maximum number of attempts
  ) {
  if (window_) {
    // clone this object nWindow times
    int t = 0;
    nWindow_ = 1;
    #ifdef MPI_H_
      MPI_Init(NULL, NULL);
      MPI_Comm_size(MPI_COMM_WORLD, &nWindow_);
      MPI_Comm_rank(MPI_COMM_WORLD, &t);
      vector<shared_ptr<WLTMMC> > clones(nWindow_);
    #endif  // MPI_H_
    #ifdef OMP_H_
      #pragma omp parallel private(t)
      {
        if (omp_get_thread_num() == 0) nWindow_ =  omp_get_num_threads();
      }
      vector<shared_ptr<WLTMMC> > clones(nWindow_);
      #pragma omp parallel private(t)
      {
        t = omp_get_thread_num();
        if (t == 0) nWindow_ = omp_get_num_threads();
    #endif  // OMP_H_

    #ifdef MPI_H_
      MPI_Barrier(MPI_COMM_WORLD);
    #endif  // MPI_H_
    #ifdef OMP_H_
      #pragma omp barrier
    #endif  // OMP_H_
    clones[t] = this->cloneShrPtr();
    clones[t]->initWindows(0);        // turn off windowing of clones

    // append output files with processor number
    stringstream ss;
    ss << "p" << t;
    clones[t]->appendFileNames(ss.str().c_str());

    // set macrostate range, and ensure the current macrostate is within range
    const double binWidth = clones[t]->c()->mBin();
    if (nOverlap_ == -1) {
      // do nothing, all windows equivalent to original
    } else if (c_->mType().compare("pairOrder") == 0) {
      vector<vector<double> > w = nWindowGrowth
        (c_->mMin(), c_->mMax(), nExp_, nWindow_, binWidth, nOverlap_);
      clones[t]->c()->initBins(w[t][1], w[t][0],
        myRound((w[t][1]-w[t][0])/binWidth));
      clones[t]->pair()->setOrder(clones[t]->c()->bin2m(0));
    } else if (c_->mType().compare("beta") == 0) {
      vector<vector<double> > w = nWindowGrowth
        (c_->mMin(), c_->mMax(), nExp_, nWindow_, binWidth, nOverlap_);
      clones[t]->c()->initBins(w[t][1], w[t][0],
        myRound((w[t][1]-w[t][0])/binWidth));
      clones[t]->criteria()->betaset(clones[t]->c()->bin2m(0));
    } else if (c_->mType().compare("pressure") == 0) {
      vector<vector<double> > w = nWindowGrowth
        (c_->mMin(), c_->mMax(), nExp_, nWindow_, binWidth, nOverlap_);
      clones[t]->c()->initBins(w[t][1], w[t][0],
                               myRound((w[t][1]-w[t][0])/binWidth));
      clones[t]->criteria()->pressureset(clones[t]->c()->bin2m(0));
    } else if (c_->mType().compare("lnpres") == 0) {
      vector<vector<double> > w = nWindowGrowth
        (c_->mMin(), c_->mMax(), nExp_, nWindow_, binWidth, nOverlap_);
      clones[t]->c()->initBins(w[t][1], w[t][0],
                               myRound((w[t][1]-w[t][0])/binWidth));
      clones[t]->criteria()->pressureset(exp(clones[t]->c()->bin2m(0)));
    } else if (betaInc_ == 0) {
      vector<vector<int> > nMolVec = nWindow(c_->mMin()+0.5, c_->mMax()-0.5,
                                             nExp_, nWindow_, nOverlap_);
      clones[t]->c()->initBins(nMolVec[t][1]+binWidth/2.,
        nMolVec[t][0]-binWidth/2.,
        myRound((nMolVec[t][1]-nMolVec[t][0])/binWidth)+1);
    } else {
      clones[t]->c()->betaset(criteria_->beta() + t*betaInc_);
      clones[t]->c()->activset(criteria_->activ()*exp(t*lnzInc_));
    }

    // #pragma omp barrier
    #ifdef MPI_H_
      MPI_Barrier(MPI_COMM_WORLD);
    #endif  // MPI_H_
    #ifdef OMP_H_
      #pragma omp barrier
    #endif  // OMP_H_
    clones[t]->writeRestart(clones[t]->rstFileName().c_str());
    #ifdef MPI_H_
      MPI_Barrier(MPI_COMM_WORLD);
    #endif  // MPI_H_
    #ifdef OMP_H_
      #pragma omp barrier
    #endif  // OMP_H_

    // initialize confswaps
    initOverlaps(t, clones);
    #ifdef MPI_H_
      MPI_Barrier(MPI_COMM_WORLD);
    #endif  // MPI_H_
    #ifdef OMP_H_
      #pragma omp barrier
    #endif  // OMP_H_
    if ( (c_->mType().compare("pairOrder") != 0) &&
         (c_->mType().compare("pressure") != 0) &&
         (c_->mType().compare("lnpres") != 0) &&
         (c_->mType().compare("beta") != 0) ) {
      clones[t]->nMolSeekInRange();
    } else {
      if (nMolSeekTarget_ != -1) clones[t]->nMolSeek(nMolSeekTarget_, 1e8);
    }

    clones[t]->checkTrialCriteria();

    // remove configuration bias trials if nMolMax/V < thres
    if (static_cast<double>(clones[t]->c()->mMax()/
        clones[t]->space()->vol()) < densThresConfigBias_) {
      cout << "t " << t << " removing configbias because nMolMax/V "
           << clones[t]->c()->mMax()/clones[t]->space()->vol()
           << " is less then thres " << densThresConfigBias_ << endl;
      clones[t]->removeConfigBias();
    }

    // initialize all clones before running
    clones[t]->zeroStat();

    // write restart files
    if (t == 0) writeRestart(rstFileName_.c_str());

    runNumSweepsExec(t, nSweeps, clones);

    #ifdef OMP_H_
      }
    #endif  // OMP_H_
  } else {
    while ( (!c_->collect() && (wlFlatTerm_ == -1) ) ||
            ( (c_->nSweep() < nSweeps) && ( (nprMax <= 0) ||
              (nAttempts_ < nprMax) ) && (wlFlatTerm_ == -1) ) ||
            ( (c_->wlFlat() < nSweeps) && ( (nprMax <= 0) ||
              (nAttempts_ < nprMax) ) && (wlFlatTerm_ != -1) )
          ) {
      attemptTrial();
    }
  }
}

/**
 * execute clones
 */
void WLTMMC::runNumSweepsExec(const int t,    //!< thread
  const int nSweeps,  //!<
  vector<shared_ptr<WLTMMC> > &clones
  ) {
  // decide whether to use MPI or OMP
  //  MPI uses restart files to communicate positions
  #ifdef MPI_H_

    // error check
    ASSERT(nFreqRestart_ == nFreqColMat_, "mpi windows use restart files,"
           << "but nFreqRestart(" << nFreqRestart_ << ") != nFreqColMat("
           << nFreqColMat_ << ")");

    // in order to print aggregate lnpi, record vector of pointers to criteria
    // initialize number of lines in criteria restart files for error check
    vector<int> nCritLines(nWindow_);
    for (int tt = 0; tt < nWindow_; ++tt) {
      stringstream ss;
      ss << rstFileName_ << "p" << tt << "criteria";
      nCritLines[tt] = numLines(ss.str().c_str());
    }

    // run until every window has atleast nSweeps
    //  sweeps updates every nFreqColMat_, so check after every nFreqColMat_
    vector<shared_ptr<CriteriaWLTMMC> > c(nWindow_);
    bool allSwept = false;
    while (allSwept == false) {
      // attempt nFreqColMat_ trials
      for (long long i = 0; i < nFreqColMat_; ++i) {
        clones[t]->attemptTrial();
      }

      // check number of sweeps from restart files
      allSwept = true;
      for (int tt = 0; tt < nWindow_; ++tt) {
        stringstream ss;
        ss << rstFileBaseName_ << "p" << tt << "criteria";
        if (myFileExists(ss.str().c_str())) {
          if (wlFlatTerm_ == -1) {
            string strtmp = fstos("nSweep", ss.str().c_str());
            string strtmp2 = fstos("collect", ss.str().c_str());
            if (!strtmp.empty() && !strtmp2.empty()) {
              if ( (stoi(strtmp) < nSweeps) || (stoi(strtmp2) == 0) ) {
                allSwept = false;
              }
            } else {
              allSwept = false;
            }
          } else {
            string strtmp = fstos("wlFlat", ss.str().c_str());
            if (!strtmp.empty()) {
              if (stoi(strtmp) < nSweeps) {
                allSwept = false;
              }
            } else {
              allSwept = false;
            }
          }
        } else {
          allSwept = false;
        }
      }

      // if master thread, use restart objects to create all threads
      // in memory of master, in order to splice
      if (t == 0) {
        bool fail = false;
        stringstream ss;
        ss << rstFileBaseName_ << "p0criteria";
        c[0] = make_shared<CriteriaWLTMMC>(ss.str().c_str());
        for (int tt = 1; tt < nWindow_; ++tt) {
          ss.str("");
          ss << rstFileBaseName_ << "p" << tt << "criteria";
          if (myFileExists(ss.str().c_str())) {
            // copy file
            std::ifstream src(ss.str().c_str(), std::fstream::binary);
            ss << "cpy";
            { std::ofstream dst(ss.str().c_str(),
              std::fstream::trunc|std::fstream::binary);
              dst << src.rdbuf();
            }

            // check number of lines if file is complete
            const int n = numLines(ss.str().c_str());
            if (n == nCritLines[tt]) {
              c[tt] = make_shared<CriteriaWLTMMC>(ss.str().c_str());
            } else {
              fail = true;
            }
          }
        }

        // print aggregate collection matrix
        if (!fail) {
          c_->spliceWindows(c);
          c_->printCollectMat(colMatFileName_.c_str());
        }
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // splice and print final colMat
    if (t == 0) {
      for (int tt = 1; tt < nWindow_; ++tt) {
        stringstream ss;
        ss << rstFileBaseName_ << "p" << tt << "criteria";
        if (myFileExists(ss.str().c_str())) {
          c[tt] = make_shared<CriteriaWLTMMC>(ss.str().c_str());
        }
      }
      c_->spliceWindows(c);
      c_->printCollectMat(colMatFileName_.c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
  #endif  // MPI_H_
  #ifdef OMP_H_

    vector<CriteriaWLTMMC*> cloneCrit(nWindow_);
    bool allSwept = false;
    for (int tt = 0; tt < nWindow_; ++tt) {
      cloneCrit[tt] = clones[tt]->c();
    }
    // #pragma omp barrier

    while (allSwept == false) {
      // attempt nFreqColMat_ trials
      for (long long i = 0; i < nFreqColMat_; ++i) {
        clones[t]->attemptTrial();
      }

      // terminate if all clones have atleast nSweeps
      if (wlFlatTerm_ == -1) {
        if (c_->minNSweep(cloneCrit) >= nSweeps) allSwept = true;
      } else {
        if (c_->minNwlFlat(cloneCrit) >= nSweeps) allSwept = true;
      }

      // print aggregate collection matrix
      if (t == 0) {
        c_->spliceWindows(cloneCrit);
        c_->printCollectMat(colMatFileName_.c_str());
      }
    }
    #pragma omp barrier
    if (t == 0) {
      c_->spliceWindows(cloneCrit);
      c_->printCollectMat(colMatFileName_.c_str());
    }
    #pragma omp barrier
  #endif  // OMP_H_
}

/**
 * resize the nmol window in WLTMMC criteria
 *  truncates liquid peak after lnPI drops as n increases by amout liquidDrop
 *  rounds to the nearest larger multiple of nround
 */
void WLTMMC::nMolResizeWindow(
  const double liquidDrop,  //!< targeted drop off of liquid peak
  const int round   //!< round up to the nearest factor
  ) {
  const int nmax = c_->nMolResizeWindow(liquidDrop, round);

  // if nMol > nmx, seek nmx
  if (space_->nMol() > nmax) nMolSeek(nmax, "", 1e7);

  const double mMax2 = nmax + 0.5;
  c_->initBins(mMax2, c_->mMin(), mMax2 - c_->mMin());
}

/**
 * seek particle number which is in the range of WLTMMC
 */
void WLTMMC::nMolSeekInRange(const int nMin,
                             const int nMax
  ) {
  double mMin = c_->mMin(), mMax = c_->mMax();
  if (nMin != -1) mMin = static_cast<double>(nMin) - 0.5;
  if (nMax != -1) mMax = static_cast<double>(nMax) + 0.5;
  if (space_->nMol() > mMax) {
    nMolSeek(myRound(mMax-0.5), 1e12);
  } else if (space_->nMol() < mMin) {
    nMolSeek(myRound(mMin+0.5), 1e12);
  }
}

/**
 * check that criteria of MC and WLTMMC match
 */
int WLTMMC::checkTrialCriteria() {
  int match = 1;
  if (c_ != criteria_) {
    match = 0;
  }
  ASSERT(match, "wltmmc criteria does not match mc criteria");
  return MC::checkTrialCriteria()*match;
}

/**
 * print saturation summary
 */
vector<double> WLTMMC::printSat() {
  std::ofstream log_(logFileName_.c_str(),
                     std::ofstream::out | std::ofstream::app);
  c_->findSat();
  c_->printRWinit();
  c_->printCollectMat(colMatFileName_.c_str());
  vector<CriteriaWLTMMC> cvec = c_->phaseSplit(c_->lnPIrw());
  const double pv = (-c_->lnPIrw().front() + log(cvec[0].lnPIarea()))
                    /space_->vol()/criteria_->beta();
  const double pl = (-c_->lnPIrw().front() + log(cvec[1].lnPIarea()))
                    /space_->vol()/criteria_->beta();
  cvec[0].lnPInorm();
  cvec[1].lnPInorm();
  log_ << "#" << className_ << " coexistance between vapor(p=" << pv
       << ", rho=" << cvec[0].lnPIaverage()/space_->vol()
       << ") and liquid(p=" << pl << ", rho="
       << cvec[1].lnPIaverage()/space_->vol() << ") finalized at z("
       << log(c_->activ()) << ", rw=" << log(c_->activrw()) << ")" << endl;
  cout << "# rhov rhol psat lnzsat phaseb" << endl;
  cout << cvec[0].lnPIaverage()/space_->vol() << " "
       << cvec[1].lnPIaverage()/space_->vol() << " " << pv << " "
       << log(c_->activrw()) << " " << cvec[0].lastbin2m() << endl;
  vector<double> returnVec;
  returnVec.push_back(cvec[0].lnPIaverage()/space_->vol());
  returnVec.push_back(cvec[1].lnPIaverage()/space_->vol());
  returnVec.push_back(pv);
  returnVec.push_back(log(c_->activrw()));
  return returnVec;
}

/**
 * write restart file
 */
void WLTMMC::writeRestart(const char* fileName) {
  MC::writeRestart(fileName);
  std::ofstream file(fileName, std::ofstream::out | std::ofstream::app);

  file << "# colMatFileName " << colMatFileName_ << endl;
  file << "# nFreqColMat " << nFreqColMat_ << endl;
  if (window_) {
    file << "# nWindow " << nWindow_ << endl;
    file << "# nExp " << nExp_ << endl;
    file << "# nOverlap " << nOverlap_ << endl;
    file << "# betaInc " << betaInc_ << endl;
    file << "# lnzInc " << lnzInc_ << endl;
  }
  file << "# densThresConfigBias " << densThresConfigBias_ << endl;

  // write restart for radial distribution functions
  if ( (nFreqGR_ != 0) && (gr_.size() > 0) ) {
    file << "# nFreqGR " << nFreqGR_ << endl;
    file << "# GRFileName " << GRFileName_ << endl;
    file << "# numGRFiles " << gr_.size() << endl;
    for (unsigned int itype = 0; itype < gr_.size(); ++itype) {
      for (unsigned int n = 0; n < gr_[itype].size(); ++n) {
        stringstream ss;
        ss << fileName << GRFileName_ << "i" << itype << "n" << n;
        if (gr_[itype][n] != NULL) gr_[itype][n]->writeRestart
          (ss.str().c_str());
      }
    }
  }

  if (wlFlatProd_ != -1) file << "# wlFlatProd " << wlFlatProd_ << endl;
  if (wlFlatTerm_ != -1) file << "# wlFlatTerm " << wlFlatTerm_ << endl;
}

/**
 * restart and run for a number of sweeps
 */
void WLTMMC::runNumSweepsRestart(
  const int nSweeps,    //!< target number of sweeps
  const char* fileName  //!< restart file name
  ) {
  // clone this object nWindow times
  int t = 0;
  nWindow_ = 1;
  #ifdef MPI_H_
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &nWindow_);
    MPI_Comm_rank(MPI_COMM_WORLD, &t);
    vector<shared_ptr<WLTMMC> > clones(nWindow_);
    MPI_Barrier(MPI_COMM_WORLD);
  #endif  // MPI_H_
  #ifdef OMP_H_
    #pragma omp parallel private(t)
    {
      if (omp_get_thread_num() == 0) nWindow_ =  omp_get_num_threads();
    }
    vector<shared_ptr<WLTMMC> > clones(nWindow_);
    #pragma omp parallel private(t)
    {
      t = omp_get_thread_num();
      nWindow_ = omp_get_num_threads();
      #pragma omp barrier
  #endif  // OMP_H_

  // read restart files
  stringstream ss;
  ss << fileName << "p" << t;
  clones[t] = make_shared<WLTMMC>(ss.str().c_str());
  clones[t]->writeRestart(clones[t]->rstFileName().c_str());

  #ifdef MPI_H_
    MPI_Barrier(MPI_COMM_WORLD);
  #endif  // MPI_H_
  #ifdef OMP_H_
    #pragma omp barrier
  #endif  // OMP_H_

  // initialize confswaps
  initOverlaps(t, clones);

  // for a more perfect restart, check if Trial parameters should be tuned
  if (nFreqRestart_ % nFreqTune_ == 0) {
    clones[t]->tuneTrialParameters();
  }
  
  runNumSweepsExec(t, nSweeps, clones);

  #ifdef OMP_H_
    }
  #endif  // OMP_H_
}

/**
 * initialize overlapping processors
 */
void WLTMMC::initOverlaps(const int t,    //!< thread
  vector<shared_ptr<WLTMMC> > &clones
  ) {
  // if configuration swap trial move exists, initialize the overlapping regions
  if (clones[t]->trialConfSwapVec_.size() == 1) {
    #ifdef MPI_H_
      TrialConfSwapMPI* trial = NULL;
      trial->initProc(t);
    #endif  // MPI_H_
    #ifdef OMP_H_
      TrialConfSwapOMP* trial = NULL;
    #endif  // OMP_H_
    trial = clones[t]->trialConfSwap(0);

    // if betaInc == 0, window density for single isotherm
    if (betaInc_ == 0) {
      // manually add smallest bins as overlapping with earlier processor,
      // skip first window
      if (t != 0) {
        for (int bin = 0; bin < nOverlap_ + 1; ++bin) {
          const double order = clones[t]->c()->bin2m(bin);
          #ifdef MPI_H_
            trial->addProcOverlap(order, t - 1);
          #endif  // MPI_H_
          #ifdef OMP_H_
            trial->addProcOverlap(order, clones[t-1]->trialConfSwap(0));
          #endif  // OMP_H_
        }
      }

      // manually add largest bins as overlapping with next processor,
      // skip last window
      if (t != nWindow_ -1) {
        const int lastbin = clones[t]->c()->nBin() - 1;
        for (int bin = lastbin - nOverlap_; bin <= lastbin; ++bin) {
          const double order = clones[t]->c()->bin2m(bin);
          #ifdef MPI_H_
            trial->addProcOverlap(order, t + 1);
          #endif  // MPI_H_
          #ifdef OMP_H_
            trial->addProcOverlap(order, clones[t+1]->trialConfSwap(0));
          #endif  // OMP_H_
        }
      }

    // if betaInc_ != 0, each processor is an isotherm of same density range
    } else {
      for (int bin = 0; bin < c_->nBin(); ++bin) {
        const double order = clones[t]->c()->bin2m(bin);
        if (t != 0) {
          #ifdef MPI_H_
            trial->addProcOverlap(order, t - 1, -betaInc_, -lnzInc_);
          #endif  // MPI_H_
          #ifdef OMP_H_
            trial->addProcOverlap(order, clones[t-1]->trialConfSwap(0),
              -betaInc_, -lnzInc_);
          #endif  // OMP_H_
        }
        if (t != nWindow_ -1) {
          #ifdef MPI_H_
            trial->addProcOverlap(order, t + 1, betaInc_, lnzInc_);
          #endif  // MPI_H_
          #ifdef OMP_H_
            trial->addProcOverlap(order, clones[t+1]->trialConfSwap(0),
              betaInc_, lnzInc_);
          #endif  // OMP_H_
        }
      }
    }
  }
}


