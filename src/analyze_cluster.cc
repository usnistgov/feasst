/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./analyze_cluster.h"

namespace feasst {

AnalyzeCluster::AnalyzeCluster(Pair *pair, const argtype &args)
  : Analyze(pair, args) {
  defaultConstruction_();
  argparse_.initArgs(className_, args);
  argparse_.checkAllArgsUsed();
}

AnalyzeCluster::AnalyzeCluster(
  Pair *pair,
  const char* fileName)
    : Analyze(pair, fileName) {
  defaultConstruction_();
  clusterCut_ = fstod("clusterCut", fileName);
}

void AnalyzeCluster::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  file << "# clusterCut " << clusterCut_ << endl;
}

void AnalyzeCluster::defaultConstruction_() {
  className_.assign("AnalyzeCluster");
  verbose_ = 0;
  clusterCut_ = 0;
  zbin = 1./(100);
  initPercolation();  //!< sets percFlag_
  nematic_.resize(3, vector<AccumulatorVec>(3));
}

void AnalyzeCluster::update(const int iMacro) {
  pair_->updateClusters(clusterCut_);
  nClusterAccVec_.accumulate(iMacro, space()->nClusters());
  const vector<int> cl = space()->clusterSizes();
  const int largestCluster = *std::max_element(cl.begin(), cl.end());
  largestClusAccVec_.accumulate(iMacro, largestCluster);
//  cout << "updating " << iMacro << endl;

  // initialize orientation histograms
  if (!space()->sphereSymMol()) {
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
  nClusSize_[iMacro].accumulate(static_cast<double>(space()->nMol())
    / static_cast<double>(space()->nClusters()));

  // compute nematic order parameter (quanternion only)
  if (space()->eulerFlag() != 1) {
    for (int iMol = 0; iMol < space()->nMol(); ++iMol) {
      const int ipart = space()->mol2part()[iMol];
      vector<double> ri(space()->dimen());
      for (int dim = 0; dim < space()->dimen(); ++dim) {
        ri[dim] = space()->x(ipart + 1, dim) - space()->x(ipart, dim);
      }
      const vector<vector<double> > mat = outerProd(ri, ri);
      for (int idim = 0; idim < space()->dimen(); ++idim) {
        for (int jdim = 0; jdim < space()->dimen(); ++jdim) {
          nematic_[idim][jdim].accumulate(iMacro, mat[idim][jdim]);
        }
      }
    }
  }

//  cout << "computing average" << endl;
  // compute the average coordination number from the contact map
  // also, compute the average z-vector orientation of contacting particles
  const vector<vector<int> > &contact = pair_->contact();
  if (contact.size() > 0) {
    std::fill(iMol2coord_.begin(), iMol2coord_.end(), 0);
    iMol2coord_.resize(space()->nMol());
    double coord = 0;
    for (int iMol = 0; iMol < space()->nMol(); ++iMol) {
      for (unsigned int index = 0; index < contact[iMol].size(); ++index) {
        const int jMol = contact[iMol][index];
        ASSERT(iMol != jMol, "iMol and jMol cannot be equal in contact map");
        coord += 1;
        ++iMol2coord_[iMol];

        // compute the z-vector orientation.
        //  place a vector along the z-axis of the reference particle,
        //  rotate and take dot product
        double cosa = -1;
        if (space()->eulerFlag() == 1) {
          vector<double> zref(space()->dimen(), 0.);
          zref[space()->dimen()-1] = 1.;
          vector<double> qi = space()->qMol(iMol);
          if (qi.size() == 4) qi.pop_back();
          vector<double> qj = space()->qMol(jMol);
          if (qj.size() == 4) qj.pop_back();
          vector<vector<double> > ri = Euler2RotMat(qi),
                                  rj = Euler2RotMat(qj);
          vector<double> zi = matVecMul(ri, zref),
                         zj = matVecMul(rj, zref);
          cosa = vecDotProd(zi, zj);

        // if not using euler angles, compute orientation by assuming that
        // it is a solid of revolution where the orientation is given by the
        // vector connecting the first particle in the molecule with the next
        } else if (!space()->sphereSymMol()) {
          vector<double> ri(space()->dimen()), rj = ri;
          const int ipart = space()->mol2part()[iMol];
          const int jpart = space()->mol2part()[jMol];
          for (int dim = 0; dim < space()->dimen(); ++dim) {
            ri[dim] = space()->x(ipart + 1, dim) - space()->x(ipart, dim);
            rj[dim] = space()->x(jpart + 1, dim) - space()->x(jpart, dim);
          }
          cosa = vecDotProd(ri, rj);
        }

        zOrient_[iMacro].accumulate(fabs(cosa));
        pcos_.accumulate(iMacro, fabs(cosa));
      }
    }
    coordNumAccVec_.accumulate(iMacro,
      coord/static_cast<double>(space()->nMol()) );
//    cout << " av " << coord/static_cast<double>(space()->nMol()) << " mac " << iMacro << endl;
  }
  // cout << "updated " << iMacro << endl;

  // test for percolation
  ASSERT( (percFlag_ >= 0) && (percFlag_ <= 2),
    "unrecognized percolation flag: " << percFlag_);
  if (percFlag_ == 1) {
    // replicate the system in each dimension
    // and see if the number of clusters change by 2^dim
    // if so, the clusters are not percolating
    shared_ptr<Space> spaceBig = space()->cloneShrPtr();
    spaceBig->replicate();
    Pair* pairBig = pair_->clone(spaceBig.get());
    pairBig->addPart();
    pairBig->updateClusters(clusterCut_);
    cout << "clus: " << space()->nClusters() << " " << spaceBig->nClusters()
         << endl;
    if (pow(2, space()->dimen())*space()->nClusters() == spaceBig->nClusters()) {
      percolation_.accumulate(iMacro, 0.);
    } else {
      percolation_.accumulate(iMacro, 1.);
    }
    delete pairBig;
  } else if (percFlag_ == 2) {
    ASSERT(pair_->className() == "PairTabular", "percolation implementation"
      << "assumes that the contact list is in the compact form, as only"
      << "implemented in PairTabular (current 4/28/2017)");
    percolation_.accumulate(iMacro, space()->percolation());
  }
}

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
         << nClusterAccVec_.vec(bin).std() << " "
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
         << pcos_.vec(bin).average() << " "
         << pcos_.vec(bin).blockStdev() << " ";

      // compute and print nematic order parameter
      vector<vector<double> > mat, evectors;
      mat.resize(3, vector<double>(3));
      vector<double> evalues;
      for (int idim = 0; idim < space()->dimen(); ++idim) {
        for (int jdim = 0; jdim < space()->dimen(); ++jdim) {
          mat[idim][jdim] = 3.*nematic_[idim][jdim].vec(bin).average()*0.5;
          if (idim == jdim) mat[idim][jdim] -= 0.5;
        }
      }
// HWH: ASSERT(jacobi(mat, evalues, evectors) != 1, "jacobi error");
      const double maxEigen = *std::max_element(evalues.begin(), evalues.end());
      ss << maxEigen << " ";
      ss << endl;
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

shared_ptr<AnalyzeCluster> makeAnalyzeCluster(Pair *pair, const argtype &args) {
  return make_shared<AnalyzeCluster>(pair, args);
}

}  // namespace feasst

