#include "./analyze_cluster.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

AnalyzeCluster::AnalyzeCluster(Space *space, Pair *pair)
  : Analyze(space, pair) {
  defaultConstruction_();
}
AnalyzeCluster::AnalyzeCluster(Space *space,
  Pair *pair,
  const char* fileName)
    : Analyze(space, pair, fileName) {
  defaultConstruction_();
  clusterCut_ = fstod("clusterCut", fileName);
}

/**
 * write restart file
 */
void AnalyzeCluster::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  file << "# clusterCut " << clusterCut_ << endl;
}

/**
 * default construction
 */
void AnalyzeCluster::defaultConstruction_() {
  className_.assign("AnalyzeCluster");
  verbose_ = 0;
  clusterCut_ = 0;
  zbin = 1./(100);
  initPercolation();  //!< sets percFlag_
}

/**
 */
void AnalyzeCluster::update(const int iMacro) {
  pair_->updateClusters(clusterCut_);
  nClusterAccVec_.accumulate(iMacro, space_->nClusters());
  const vector<int> cl = space_->clusterSizes();
  const int largestCluster = *std::max_element(cl.begin(), cl.end());
  largestClusAccVec_.accumulate(iMacro, largestCluster);
//  cout << "updating " << iMacro << endl;

  // initialize orientation histograms
  if (!space_->sphereSymMol()) {
    if (iMacro >= static_cast<int>(zOrient_.size())) {
      const int zsorig = zOrient_.size();
      for (int inc = 0; inc < iMacro - zsorig + 1; ++inc) {
        Histogram h(zbin);
        // h.centerZero();
        zOrient_.push_back(h);
      }
    }
  }

  // initialize cluster size histogram
  const int sorig = nClusSize_.size();
  if (iMacro >= sorig) {
    for (int inc = 0; inc < iMacro - sorig + 1; ++inc) {
      Histogram h(1.);
      nClusSize_.push_back(h);
    }
  }
  nClusSize_[iMacro].accumulate(static_cast<double>(space_->nMol())
    / static_cast<double>(space_->nClusters()));

//  cout << "computing average" << endl;
  // compute the average coordination number from the contact map
  // also, compute the average z-vector orientation of contacting particles
  const vector<vector<int> > &contact = pair_->contact();
  if (contact.size() > 0) {
    double coord = 0;
    for (int iMol = 0; iMol < space_->nMol(); ++iMol) {
      for (unsigned int index = 0; index < contact[iMol].size(); ++index) {
        const int jMol = contact[iMol][index];
        ASSERT(iMol != jMol, "iMol and jMol cannot be equal in contact map");
        coord += 1;

        // compute the z-vector orientation.
        //  place a vector along the z-axis of the reference particle,
        //  rotate and take dot product
        double cosa = -1;
        if (space_->eulerFlag() == 1) {
          vector<double> zref(space_->dimen(), 0.);
          zref[space_->dimen()-1] = 1.;
          vector<double> qi = space_->qMol(iMol);
          if (qi.size() == 4) qi.pop_back();
          vector<double> qj = space_->qMol(jMol);
          if (qj.size() == 4) qj.pop_back();
          vector<vector<double> > ri = Euler2RotMat(qi),
                                  rj = Euler2RotMat(qj);
          vector<double> zi = matVecMul(ri, zref),
                         zj = matVecMul(rj, zref);
          cosa = vecDotProd(zi, zj);

        // if not using euler angles, compute orientation by assuming that
        // it is a solid of revolution where the orientation is given by the
        // vector connecting the first particle in the molecule with the next
        } else if (!space_->sphereSymMol()) {
          vector<double> ri(space_->dimen()), rj = ri;
          const int ipart = space_->mol2part()[iMol];
          const int jpart = space_->mol2part()[jMol];
          for (int dim = 0; dim < space_->dimen(); ++dim) {
            ri[dim] = space_->x(ipart + 1, dim) - space_->x(ipart, dim);
            rj[dim] = space_->x(jpart + 1, dim) - space_->x(jpart, dim);
          }
          cosa = vecDotProd(ri, rj);
        }

        zOrient_[iMacro].accumulate(fabs(cosa));
      }
    }
    coordNumAccVec_.accumulate(iMacro,
      coord/static_cast<double>(space_->nMol()) );
  }
  // cout << "updated " << iMacro << endl;

  // test for percolation
  ASSERT( (percFlag_ >= 0) && (percFlag_ <= 2),
    "unrecognized percolation flag: " << percFlag_);
  if (percFlag_ == 1) {
    // replicate the system in each dimension
    // and see if the number of clusters change by 2^dim
    // if so, the clusters are not percolating
    shared_ptr<Space> spaceBig = space_->cloneShrPtr();
    spaceBig->replicate();
    Pair* pairBig = pair_->clone(spaceBig.get());
    pairBig->addPart();
    pairBig->updateClusters(clusterCut_);
    cout << "clus: " << space_->nClusters() << " " << spaceBig->nClusters()
         << endl;
    if (pow(2, space_->dimen())*space_->nClusters() == spaceBig->nClusters()) {
      percolation_.accumulate(iMacro, 0.);
    } else {
      percolation_.accumulate(iMacro, 1.);
    }
    delete pairBig;
  } else if (percFlag_ == 2) {
    ASSERT(pair_->className() == "PairTabular", "percolation implementation"
      << "assumes that the contact list is in the compact form, as only"
      << "implemented in PairTabular (current 4/28/2017)");
    percolation_.accumulate(iMacro, space_->percolation());
  }
}

/**
 * printer
 */
void AnalyzeCluster::write(CriteriaWLTMMC *c) {
  // initialize output
  fileBackUp(fileName_.c_str());
  std::ofstream file(fileName_.c_str());
  stringstream ss;
  ss << "# " << c->mType() << " nClusterAv nClusterSdtev nClusterBlockStdev "
     << "coordNum coordNumStd" << endl;
  if (fileName_.empty()) {
    cout << ss.str();
  } else {
    file << ss.str();
  }

  // print
  for (int bin = 0; bin < c->nBin(); ++bin) {
    ss.str("");
    ss << c->bin2m(bin) << " ";
    if (nClusterAccVec_.size() <= bin) {
      ss << "-1 -1 -1" << endl;
    } else {
      ss << nClusterAccVec_.vec(bin).average() << " "
         << nClusterAccVec_.vec(bin).stdev() << " "
         << nClusterAccVec_.vec(bin).blockStdev() << " "
         << coordNumAccVec_.vec(bin).average() << " "
         << coordNumAccVec_.vec(bin).blockStdev() << " ";
      if (static_cast<int>(zOrient_.size()) > bin) {
        // compute relative probability of |cosa|<0.5 or >0.5
        Histogram* h = &(zOrient_[bin]);
        double count = 0., tot = 0.;
        if (h->size() > 0) {
          for (int hbin = 0; hbin < h->size(); ++hbin) {
            const double hval = h->hist()[hbin];
            tot += hval;
            if (h->bin2m(hbin) <= 0.5) count += hval;
          }
          ss << count/tot << " ";
        } else {
          ss << "0 ";
        }
      }
      ss << largestClusAccVec_.vec(bin).average() << " "
         << largestClusAccVec_.vec(bin).blockStdev() << " "
         << percolation_.vec(bin).average() << " "
         << percolation_.vec(bin).blockStdev() << " "
         << endl;
    }
    if (fileName_.empty()) {
      cout << ss.str();
    } else {
      file << ss.str();
    }
  }

  // print zOrient probability distributions, 1 file per macrostate
  if (zOrient_.size() > 0) {
    for (unsigned int iMacro = 0; iMacro < zOrient_.size(); ++iMacro) {

      // initialize output
      stringstream ssfn;
      ssfn << fileName_ << "i" << iMacro;
      fileBackUp(ssfn.str().c_str());
      std::ofstream file(ssfn.str().c_str());
      stringstream ss;
      ss << "# m p(m)" << endl;
      if (fileName_.empty()) {
        cout << ss.str();
      } else {
        file << ss.str();
      }

      // print histogram
      ss.str("");
      Histogram* h = &(zOrient_[iMacro]);
      for (int bin = 0; bin < h->size(); ++bin) {
        ss << h->bin2m(bin) << " " << h->hist()[bin] << endl;
      }
      if (fileName_.empty()) {
        cout << ss.str();
      } else {
        file << ss.str();
      }

      // initialize output
      ssfn.str("");
      ssfn << fileName_ << "histi" << iMacro;
      fileBackUp(ssfn.str().c_str());
      std::ofstream file2(ssfn.str().c_str());
      ss.str("");
      ss << "# m p(m)" << endl;
      if (fileName_.empty()) {
        cout << ss.str();
      } else {
        file2 << ss.str();
      }

      // print histogram
      ss.str("");
      h = &(nClusSize_[iMacro]);
      for (int bin = 0; bin < h->size(); ++bin) {
        ss << h->bin2m(bin) << " " << h->hist()[bin] << endl;
      }
      if (fileName_.empty()) {
        cout << ss.str();
      } else {
        file2 << ss.str();
      }
    }
  }
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

