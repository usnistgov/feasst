/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./mc_wltmmc.h"

namespace feasst {

WLTMMC::WLTMMC(Space* space,
  Pair* pair,
  CriteriaWLTMMC* criteria
  )
    : MC(space, pair, criteria),
      c_(criteria) {
  defaultConstruction_();
}
WLTMMC::WLTMMC(const char* fileName)
    : MC(fileName) {
  defaultConstruction_();

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

  strtmp = fstos("procFileAppend", fileName);
  if (!strtmp.empty()) {
    procFileAppend_ = strtmp;
  }

  strtmp = fstos("wlFlatProd", fileName);
  if (!strtmp.empty()) {
    wlFlatProd_ = stoi(strtmp);
  }
  strtmp = fstos("wlFlatTerm", fileName);
  if (!strtmp.empty()) {
    wlFlatTerm_ = stoi(strtmp);
  }

  if (c_->nSweep() > 1) production_ = 1;
}

void WLTMMC::defaultConstruction_() {
  className_.assign("WLTMMC");
  verbose_ = 0;
  nFreqColMat_ = 1e6;
  colMatFileName_.assign("colMat");
  initWindows(0);
  betaInc_ = 0;
  lnzInc_ = 0;
  densThresConfigBias_ = 0;
  nMolSeekTarget_ = -1;
  wlFlatProd_ = -1;
  wlFlatTerm_ = -1;
  setProcessorFileDescription();
}

WLTMMC::~WLTMMC() {
  // if (pair_ != NULL) pair_->checkEnergy(1e-9, 0);
  checkTrialCriteria();
  if (criteriaOwned_) {
    delete c_;
  }
}

void WLTMMC::reconstruct() {
  CriteriaWLTMMC* c = c_->clone();
  criteriaOwned_ = true;
  c_ = c;  // swap temporaries for exception safety
  criteria_ = reinterpret_cast<Criteria*>(c_);
  MC::reconstruct();
}

WLTMMC* WLTMMC::clone() const {
  WLTMMC* mc = new WLTMMC(*this);
  mc->reconstruct();
  return mc;
}

shared_ptr<MC> WLTMMC::cloneImpl() const {
  shared_ptr<WLTMMC> mc = make_shared<WLTMMC>(*this);
  mc->reconstruct();
  return mc;
}

void WLTMMC::afterAttempt_() {
  if ( ( (c_->nSweep() > 1) ||
         ((wlFlatProd_ != -1) && (c_->wlFlat() >= wlFlatProd_) ) )
    && (production_ == 0) ) {
    initProduction();
  }
  afterAttemptBase_();
  if (nAttempts_ % nFreqColMat_ == 0) {
    c_->lnPIupdate();
    if (!colMatFileName_.empty()) {
//      c_->lnPIpressureIso(space_->volume());
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

  // run analyzers
  for (vector<shared_ptr<Analyze>>::iterator it = analyzeVec_.begin();
       it != analyzeVec_.end(); ++it) {
    if (nAttempts_ % (*it)->nFreq() == 0) {
      (*it)->update(c_->iMacro());
    }
    if (nAttempts_ % (*it)->nFreqPrint() == 0) {
      (*it)->write(c_);
    }
  }
}

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
    #endif  // MPI_H_
    #ifdef _OPENMP
      #pragma omp parallel private(t)
      {
        if (omp_get_thread_num() == 0) nWindow_ =  omp_get_num_threads();
      }
    #endif  // _OPENMP
    vector<shared_ptr<WLTMMC> > clones(nWindow_);
    #ifdef _OPENMP
      #pragma omp parallel private(t)
      {
        t = omp_get_thread_num();
        if (t == 0) nWindow_ = omp_get_num_threads();
    #endif  // _OPENMP

    #ifdef MPI_H_
      MPI_Barrier(MPI_COMM_WORLD);
    #endif  // MPI_H_
    #ifdef _OPENMP
      #pragma omp barrier
    #endif  // _OPENMP
    clones[t] = this->cloneShrPtr();
    clones[t]->initWindows(0);        // turn off windowing of clones

    // append output files with processor number
    stringstream ss;
    ss << procFileAppend_ << t;
    clones[t]->appendFileNames(ss.str().c_str());

    // set macrostate range, and ensure the current macrostate is within range
    const double binWidth = clones[t]->c()->mBin();
    if (nOverlap_ == -1) {
      // do nothing, all windows equivalent to original
    } else if (c_->mType().compare("pairOrder") == 0) {
      vector<vector<double> > w = nWindowGrowth
        (c_->mMin(), c_->mMax(), nExp_, nWindow_, binWidth, nOverlap_);
      clones[t]->c()->initBins(w[t][1], w[t][0],
        feasstRound((w[t][1]-w[t][0])/binWidth));
      clones[t]->pair()->setOrder(clones[t]->c()->bin2m(0));
    } else if (c_->mType().compare("beta") == 0) {
      vector<vector<double> > w = nWindowGrowth
        (c_->mMin(), c_->mMax(), nExp_, nWindow_, binWidth, nOverlap_);
      clones[t]->c()->initBins(w[t][1], w[t][0],
        feasstRound((w[t][1]-w[t][0])/binWidth));
      clones[t]->criteria()->betaset(clones[t]->c()->bin2m(0));
    } else if (c_->mType().compare("pressure") == 0) {
      vector<vector<double> > w = nWindowGrowth
        (c_->mMin(), c_->mMax(), nExp_, nWindow_, binWidth, nOverlap_);
      clones[t]->c()->initBins(w[t][1], w[t][0],
                               feasstRound((w[t][1]-w[t][0])/binWidth));
      clones[t]->criteria()->pressureset(clones[t]->c()->bin2m(0));
    } else if (c_->mType().compare("lnpres") == 0) {
      vector<vector<double> > w = nWindowGrowth
        (c_->mMin(), c_->mMax(), nExp_, nWindow_, binWidth, nOverlap_);
      clones[t]->c()->initBins(w[t][1], w[t][0],
                               feasstRound((w[t][1]-w[t][0])/binWidth));
      clones[t]->criteria()->pressureset(exp(clones[t]->c()->bin2m(0)));
    } else if (betaInc_ == 0) {
      vector<vector<int> > nMolVec = nWindow(c_->mMin()+0.5, c_->mMax()-0.5,
                                             nExp_, nWindow_, nOverlap_);
      clones[t]->c()->initBins(nMolVec[t][1]+binWidth/2.,
        nMolVec[t][0]-binWidth/2.,
        feasstRound((nMolVec[t][1]-nMolVec[t][0])/binWidth)+1);
    } else {
      clones[t]->c()->betaset(criteria_->beta() + t*betaInc_);
      clones[t]->c()->activset(criteria_->activ()*exp(t*lnzInc_));
    }

    // #pragma omp barrier
    #ifdef MPI_H_
      MPI_Barrier(MPI_COMM_WORLD);
    #endif  // MPI_H_
    #ifdef _OPENMP
      #pragma omp barrier
    #endif  // _OPENMP
    clones[t]->writeRestart(clones[t]->rstFileName().c_str());
    #ifdef MPI_H_
      MPI_Barrier(MPI_COMM_WORLD);
    #endif  // MPI_H_
    #ifdef _OPENMP
      #pragma omp barrier
    #endif  // _OPENMP

    // initialize confswaps
    initOverlaps(t, &clones);
    #ifdef MPI_H_
      MPI_Barrier(MPI_COMM_WORLD);
    #endif  // MPI_H_
    #ifdef _OPENMP
      #pragma omp barrier
    #endif  // _OPENMP
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
        clones[t]->space()->volume()) < densThresConfigBias_) {
      cout << "t " << t << " removing configbias because nMolMax/V "
           << clones[t]->c()->mMax()/clones[t]->space()->volume()
           << " is less then thres " << densThresConfigBias_ << endl;
      clones[t]->removeConfigBias();
    }

    // initialize all clones before running
    clones[t]->zeroStat();

    // write restart files
    if (t == 0) writeRestart(rstFileName_.c_str());

    runNumSweepsExec_(t, nSweeps, &clones);

    #ifdef _OPENMP
      }
    #endif  // _OPENMP
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

void WLTMMC::runNumSweepsExec_(const int t,    //!< thread
  const int nSweeps,  //!<
  vector<shared_ptr<WLTMMC> > *clones
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
        (*clones)[t]->attemptTrial();
      }

      // check number of sweeps from restart files
      allSwept = true;
      for (int tt = 0; tt < nWindow_; ++tt) {
        stringstream ss;
        ss << rstFileBaseName_ << "p" << tt << "criteria";
        if (fileExists(ss.str().c_str())) {
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
          if (fileExists(ss.str().c_str())) {
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
        if (fileExists(ss.str().c_str())) {
          c[tt] = make_shared<CriteriaWLTMMC>(ss.str().c_str());
        }
      }
      c_->spliceWindows(c);
      c_->printCollectMat(colMatFileName_.c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
  #endif  // MPI_H_
  #ifdef _OPENMP

    vector<CriteriaWLTMMC*> cloneCrit(nWindow_);
    bool allSwept = false;
    for (int tt = 0; tt < nWindow_; ++tt) {
      cloneCrit[tt] = (*clones)[tt]->c();
    }
    // #pragma omp barrier

    while (allSwept == false) {
      // attempt nFreqColMat_ trials
      for (long long i = 0; i < nFreqColMat_; ++i) {
        (*clones)[t]->attemptTrial();
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
  #endif  // _OPENMP
}

void WLTMMC::nMolSeekInRange(const int nMin,
                             const int nMax
  ) {
  double mMin = c_->mMin(), mMax = c_->mMax();
  if (nMin != -1) mMin = static_cast<double>(nMin) - 0.5;
  if (nMax != -1) mMax = static_cast<double>(nMax) + 0.5;
  // if order parameter is number of first molecules
  if (c_->mType().compare("nmol0") == 0) {
    const int nMolType0 = space_->nMolType()[0];
    if (nMolType0 > mMax) {
      nMolSeek(space_->nMol() - nMolType0 + feasstRound(mMax - 0.5), 1e12);
    } else if (nMolType0 < mMin) {
      nMolSeek(space_->nMol() - nMolType0 + feasstRound(mMin + 0.5), 1e12);
    }
  // if order parameter is total number of molecules
  } else if (c_->mType().compare("nmol") == 0) {
    if (space_->nMol() > mMax) {
      nMolSeek(feasstRound(mMax-0.5), 1e12);
    } else if (space_->nMol() < mMin) {
      nMolSeek(feasstRound(mMin+0.5), 1e12);
    }
  } else {
    ASSERT(0, "unrecognized criteria order type");
  }
}

int WLTMMC::checkTrialCriteria() {
  int match = 1;
  if (c_ != criteria_) {
    match = 0;
  }
  ASSERT(match, "wltmmc criteria does not match mc criteria");
  return MC::checkTrialCriteria()*match;
}

vector<double> WLTMMC::printSat() {
  std::ofstream log_(logFileName_.c_str(),
                     std::ofstream::out | std::ofstream::app);
  c_->findSat();
  c_->printRWinit();
  c_->printCollectMat(colMatFileName_.c_str());
  vector<CriteriaWLTMMC> cvec = c_->phaseSplit(c_->lnPIrw());
  const double pv = (-c_->lnPIrw().front() + log(cvec[0].lnPIarea()))
                    /space_->volume()/criteria_->beta();
  const double pl = (-c_->lnPIrw().front() + log(cvec[1].lnPIarea()))
                    /space_->volume()/criteria_->beta();
  cvec[0].lnPInorm();
  cvec[1].lnPInorm();
  log_ << "#" << className_ << " coexistance between vapor(p=" << pv
       << ", rho=" << cvec[0].lnPIaverage()/space_->volume()
       << ") and liquid(p=" << pl << ", rho="
       << cvec[1].lnPIaverage()/space_->volume() << ") finalized at z("
       << log(c_->activ()) << ", rw=" << log(c_->activrw()) << ")" << endl;
  cout << "# rhov rhol psat lnzsat phaseb" << endl;
  cout << cvec[0].lnPIaverage()/space_->volume() << " "
       << cvec[1].lnPIaverage()/space_->volume() << " " << pv << " "
       << log(c_->activrw()) << " " << cvec[0].lastbin2m() << endl;
  vector<double> returnVec;
  returnVec.push_back(cvec[0].lnPIaverage()/space_->volume());
  returnVec.push_back(cvec[1].lnPIaverage()/space_->volume());
  returnVec.push_back(pv);
  returnVec.push_back(log(c_->activrw()));
  return returnVec;
}

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
  file << "# procFileAppend " << procFileAppend_ << endl;

  if (wlFlatProd_ != -1) file << "# wlFlatProd " << wlFlatProd_ << endl;
  if (wlFlatTerm_ != -1) file << "# wlFlatTerm " << wlFlatTerm_ << endl;
}

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
    MPI_Barrier(MPI_COMM_WORLD);
  #endif  // MPI_H_
  #ifdef _OPENMP
    #pragma omp parallel private(t)
    {
      if (omp_get_thread_num() == 0) nWindow_ =  omp_get_num_threads();
    }
  #endif  // _OPENMP

  vector<shared_ptr<WLTMMC> > clones(nWindow_);
  #ifdef _OPENMP
    #pragma omp parallel private(t)
    {
      t = omp_get_thread_num();
      nWindow_ = omp_get_num_threads();
      #pragma omp barrier
  #endif  // _OPENMP

  // read restart files
  stringstream ss;
  ss << fileName << "p" << t;
  clones[t] = make_shared<WLTMMC>(ss.str().c_str());
  clones[t]->writeRestart(clones[t]->rstFileName().c_str());

  #ifdef MPI_H_
    MPI_Barrier(MPI_COMM_WORLD);
  #endif  // MPI_H_
  #ifdef _OPENMP
    #pragma omp barrier
  #endif  // _OPENMP

  // initialize confswaps
  initOverlaps(t, &clones);

  // for a more perfect restart, check if Trial parameters should be tuned
  if (nFreqTune_ != 0) {
    if (nFreqRestart_ % nFreqTune_ == 0) {
      clones[t]->tuneTrialParameters_();
    }
  }

  // monkey patch
  for (vector<shared_ptr<Analyze>>::iterator it = analyzeVec_.begin();
       it != analyzeVec_.end(); ++it) {
    (*it)->modifyRestart(clones[t]);
  }

  runNumSweepsExec_(t, nSweeps, &clones);

  #ifdef _OPENMP
    }
  #endif  // _OPENMP
}

void WLTMMC::initOverlaps(const int t,    // thread
  vector<shared_ptr<WLTMMC> > *clones
  ) {
  // if configuration swap trial move exists, initialize the overlapping regions
  #if defined (MPI_H_) || (_OPENMP)
  if ((*clones)[t]->trialConfSwapVec_.size() == 1) {
    #ifdef MPI_H_
      TrialConfSwapMPI* trial = NULL;
      trial->initProc(t);
    #endif  // MPI_H_
    #ifdef _OPENMP
      TrialConfSwapOMP* trial = NULL;
    #endif  // _OPENMP
    trial = (*clones)[t]->trialConfSwap(0);

    // if betaInc == 0, window density for single isotherm
    if (betaInc_ == 0) {
      // manually add smallest bins as overlapping with earlier processor,
      // skip first window
      if (t != 0) {
        for (int bin = 0; bin < nOverlap_ + 1; ++bin) {
          const double order = (*clones)[t]->c()->bin2m(bin);
          #ifdef MPI_H_
            trial->addProcOverlap(order, t - 1);
          #endif  // MPI_H_
          #ifdef _OPENMP
            trial->addProcOverlap(order, (*clones)[t-1]->trialConfSwap(0));
          #endif  // _OPENMP
        }
      }

      // manually add largest bins as overlapping with next processor,
      // skip last window
      if (t != nWindow_ -1) {
        const int lastbin = (*clones)[t]->c()->nBin() - 1;
        for (int bin = lastbin - nOverlap_; bin <= lastbin; ++bin) {
          const double order = (*clones)[t]->c()->bin2m(bin);
          #ifdef MPI_H_
            trial->addProcOverlap(order, t + 1);
          #endif  // MPI_H_
          #ifdef _OPENMP
            trial->addProcOverlap(order, (*clones)[t+1]->trialConfSwap(0));
          #endif  // _OPENMP
        }
      }

    // if betaInc_ != 0, each processor is an isotherm of same density range
    } else {
      for (int bin = 0; bin < c_->nBin(); ++bin) {
        const double order = (*clones)[t]->c()->bin2m(bin);
        if (t != 0) {
          #ifdef MPI_H_
            trial->addProcOverlap(order, t - 1, -betaInc_, -lnzInc_);
          #endif  // MPI_H_
          #ifdef _OPENMP
            trial->addProcOverlap(order, (*clones)[t-1]->trialConfSwap(0),
              -betaInc_, -lnzInc_);
          #endif  // _OPENMP
        }
        if (t != nWindow_ -1) {
          #ifdef MPI_H_
            trial->addProcOverlap(order, t + 1, betaInc_, lnzInc_);
          #endif  // MPI_H_
          #ifdef _OPENMP
            trial->addProcOverlap(order, (*clones)[t+1]->trialConfSwap(0),
              betaInc_, lnzInc_);
          #endif  // _OPENMP
        }
      }
    }
  }
  #endif  // MPI_H_ || _OPENMP
}

}  // namespace feasst
