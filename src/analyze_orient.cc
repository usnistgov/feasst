#include "./analyze_orient.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

AnalyzeOrient::AnalyzeOrient(Space *space, Pair *pair)
  : Analyze(space, pair) {
  defaultConstruction();
}
AnalyzeOrient::AnalyzeOrient(Space *space,
  Pair *pair,
  const char* fileName)
  : Analyze(space, pair, fileName) {
  defaultConstruction();
  // clusterCut_ = fstod("clusterCut", fileName);
}

/**
 * write restart file
 */
void AnalyzeOrient::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  // file << "# clusterCut " << clusterCut_ << endl;
}

/**
 * default construction
 */
void AnalyzeOrient::defaultConstruction() {
  className_.assign("AnalyzeOrient");
  verbose_ = 0;
  zbin = 0.05;
  ASSERT(!space_->sphereSymMol(),
    "cannot compute orientation properties of spherically symmetric particles");
  ASSERT(space_->eulerFlag(), "implemented only for euler angles");
}

/**
 */
void AnalyzeOrient::update(const int iMacro) {
  // cout << "updating " << iMacro << endl;
  // initialize histogram
  if (h_.size() <= 0) {
    h_.initBinWidth(zbin);
    h_.accumulate(zbin/2.);
    h_.accumulate(space_->minl());
  }

  // initialize orientation accumulator
  if (iMacro >= static_cast<int>(zOrient_.size())) {
    const int zsorig = zOrient_.size();
    for (int inc = 0; inc < iMacro - zsorig + 1; ++inc) {
      AccumulatorVec accVec;
      zOrient_.push_back(accVec);
    }
  }

  // loop through each pair of molecules,
  // compute center separation distance, and angle between the axis of symmetry
  const int nMol = space_->nMol();
  const vector<double> &x = space_->x();
  const vector<int> &mol2part = space_->mol2part();
  const int dimen_ = space_->dimen();
  for (int iMol = 0; iMol < nMol - 1; ++iMol) {
    const int ipart = mol2part[iMol];
    for (int jMol = iMol + 1; jMol < nMol; ++jMol) {
      const int jpart = mol2part[jMol];

      // separation vector, xij with periodic boundary conditions
      vector<double> xij(dimen_);
      for (int i = 0; i < dimen_; ++i) {
        xij[i] = x[dimen_*ipart+i] - x[dimen_*jpart+i];
      }
      const vector<double> dx = space_->pbc(xij);
      for (int dim = 0; dim < dimen_; ++dim) {
        xij[dim] += dx[dim];
      }
      const double r2 = vecDotProd(xij, xij);
      const double r = sqrt(r2);

      // compute the orientation
      vector<double> iv(dimen_), jv(dimen_);
      for (int dim = 0; dim < dimen_; ++dim) {
        iv[dim] = x[dimen_*(ipart+1)+dim] - x[dimen_*ipart+dim];
        jv[dim] = x[dimen_*(jpart+1)+dim] - x[dimen_*jpart+dim];
      }
      const double cosa = vecDotProd(iv, jv);
//      const double ti = vecDotProd(iv, iv),
//                   tj = vecDotProd(jv, jv);
//      ASSERT(fabs(ti-1)<1e-10, "testi failed");
//      ASSERT(fabs(tj-1)<1e-10, "testj failed");
      // cout << "r " << r << " bin " << h_.bin(r) << " cosa " << cosa << endl;
      zOrient_[iMacro].accumulate(h_.bin(r), fabs(cosa));
    }
  }
}

/**
 * printer
 */
void AnalyzeOrient::write() {
  // print zOrient probability distributions, 1 file per macrostate
  if (zOrient_.size() > 0) {
    for (unsigned int iMacro = 0; iMacro < zOrient_.size(); ++iMacro) {
      // initialize output
      stringstream ssfn;
      ssfn << fileName_ << "i" << iMacro;
      fileBackUp(ssfn.str().c_str());
      std::ofstream file(ssfn.str().c_str());
      stringstream ss;
      ss << "#r |cos|" << endl;
      if (fileName_.empty()) {
        cout << ss.str();
      } else {
        file << ss.str();
      }

      // print distribution
      ss.str("");
      AccumulatorVec* av = &(zOrient_[iMacro]);
      for (int bin = 0; bin < h_.size(); ++bin) {
        if (bin < av->size()) {
          ss << h_.bin2m(bin) << " " << (av->vec(bin)).average() << endl;
          // if (av->size() zOrien<< h_->hist()[bin] << endl;
        }
      }
      if (fileName_.empty()) {
        cout << ss.str();
      } else {
        file << ss.str();
      }
    }
  }
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_


