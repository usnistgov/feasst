#include <limits>
#include <algorithm>
#include "./space.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

Space::Space(int dimen, int id)
  : dimen_(dimen),
    id_(id) {
  className_.assign("Space");
  defaultConstruction_();
}
Space::Space(const char* fileName) {
  className_.assign("Space");
  ASSERT(fileExists(fileName),
         "restart file(" << fileName << ") doesn't exist");

  // cout << " initialize identity and dimensionality, plus defaults" << endl;
  id_ = fstoi("id", fileName);
  dimen_ = fstoi("dimen", fileName);
  defaultConstruction_();

  // cout << " initialize simulation domain" << endl;
  for (int dim = 0; dim < dimen_; ++dim) {
    stringstream ss;
    ss << "l" << dim;
    l_[dim] = fstod(ss.str().c_str(), fileName);
  }
  string strtmp = fstos("xyTilt", fileName);
  if (!strtmp.empty()) xyTilt_ = stod(strtmp);
  strtmp = fstos("xzTilt", fileName);
  if (!strtmp.empty()) xzTilt_ = stod(strtmp);
  strtmp = fstos("yzTilt", fileName);
  if (!strtmp.empty()) yzTilt_ = stod(strtmp);
  strtmp = fstos("floppyBoxFlag", fileName);
  if (!strtmp.empty()) floppyBox_ = stoi(strtmp);
  strtmp = fstos("maxlFlag", fileName);
  if (!strtmp.empty()) {
    maxlFlag_ = stoi(strtmp);
    maxl_.resize(dimen_);
    for (int dim = 0; dim < dimen_; ++dim) {
      stringstream ss;
      ss << "maxboxl" << dim;
      maxl_[dim] = fstod(ss.str().c_str(), fileName);
    }
  }

  // addmolinits
  int nMolTypes = fstoi("naddmolinits", fileName);
  for (int i = 0; i < nMolTypes; ++i) {
    stringstream ss;
    ss << "addmolinittype" << i;
    const string addmolstr = fstos(ss.str().c_str(), fileName);
    addMolInit(addmolstr.c_str());
  }

  // cout << " initialize molecular descriptions" << endl;
  nMolTypes = fstoi("nMolTypes", fileName);
  for (int i = 0; i < nMolTypes; ++i) {
    stringstream ss;
    ss << "moltype" << i;
    const string addmolstr = fstos(ss.str().c_str(), fileName);
    ss.str("");
    ss << "nmoloft" << i;
    const int naddmol = fstoi(ss.str().c_str(), fileName);
    // addMolInit(addmolstr.c_str());
    for (int iMol = 0; iMol < naddmol; ++iMol) {
      addMol(addmolstr.c_str());
    }
  }

  // cout << " update tags" << endl;
  strtmp = fstos("ntags", fileName);
  if (!strtmp.empty()) {
    const int ntags = fstoi("ntags", fileName);
    for (int i = 0; i < ntags; ++i) {
      stringstream ss;
      ss << "tag" << i;
      tagAtom(fstoi(ss.str().c_str(), fileName));
    }
  }

  // print something about cluster analysis
  const int nClusterTypes = fstoi("nClusterTypes", fileName);
  for (int i = 0; i < nClusterTypes; ++i) {
    stringstream ss;
    ss << "clusterType" << i;
    const int type = fstoi(ss.str().c_str(), fileName);
    clusterType_.push_back(type);
  }

  // determine whether to use euler angles or quaternions
  strtmp = fstos("eulerFlag", fileName);
  if (!strtmp.empty()) {
    eulerFlag_ = stoi(strtmp);
  }

  strtmp = fstos("equiMolarFlag", fileName);
  if (!strtmp.empty()) {
    equiMolar_ = stoi(strtmp);
  }

  // cout << " open file and skip header lines" << endl;
  std::ifstream fs(fileName);
  string line;
  const int nLines = numLines(fileName);
  int nSkip;
  if (sphereSymMol_) {
    nSkip = nLines - natom();
  } else {
    nSkip = nLines - (natom() + nMol());
  }
  for (int i = 0; i < nSkip; ++i) getline(fs, line);

  // read configuration
  if (sphereSymMol_) {
    // cout << " read particle positions" << endl;
    for (int i = 0; i < natom(); ++i) {
      for (int dim = 0; dim < dimen_; ++dim) {
        fs >> x_[dimen_*i+dim];
      }
      getline(fs, line);
    }
  } else {
    // for each molecule, read position of first(pivot) atom, xMolRef, and qMol
    for (int iMol = 0; iMol < nMol(); ++iMol) {
      const int iAtom = mol2part_[iMol];
      for (int dim = 0; dim < dimen_; ++dim) {
        fs >> x_[dimen_*iAtom + dim];
      }
      getline(fs, line);
      for (int dim = 0; dim < dimen_; ++dim) {
        xMolRef_[iMol][0][dim] = 0.;
      }
      for (unsigned int ipart = 1; ipart < xMolRef_[iMol].size(); ++ipart) {
        for (int dim = 0; dim < dimen_; ++dim) {
          fs >> xMolRef_[iMol][ipart][dim];
        }
        getline(fs, line);
      }
      for (int qdim = 0; qdim < qdim_; ++qdim) {
        fs >> qMol_[qdim_*iMol + qdim];
      }
      getline(fs, line);
      quat2pos(iMol);
    }
  }

//  // cout << " initialize quaternions and xMolRef" << endl;
//  if ( (!sphereSymMol_) && (natom() != 1) ) qMolInit();

  // cout << " initialize cell list" << endl;
  cellAtomCut_ = fstoi("cellAtomCut", fileName);
  initCellAtomCut(cellAtomCut_);
  strtmp = fstos("cellType", fileName);
  if (!strtmp.empty()) {
    cellType_ = stoi(strtmp);
    dCellMin_ = fstod("dCellMin", fileName);
    updateCells(dCellMin_);
  } else {
    cellType_ = 0;
  }

  // initialize random number generator
  initRNG(fileName);
}

void Space::defaultConstruction_() {
  verbose_ = 0;
  fastDel_ = false;
  fastDelMol_ = -1;
  if (dimen_ == 3) {
    qdim_ = dimen_ + 1;
  } else if (dimen_ == 2) {
    qdim_ = 1;
  } else {
    ASSERT(dimen_ == 1, "this dimensionality(" << dimen_ << ") is not"
           << "supported");
  }
  l_.resize(dimen_);
  xyTilt_ = 0.;
  xzTilt_ = 0.;
  yzTilt_ = 0.;
  floppyBox_ = 0;
  maxlFlag_ = 0;
  mol2part_.push_back(0);
  nType_.resize(1, 0);
  nMolType_.resize(1, 0);
  sphereSymMol_ = true;
  cellOff();
  initCellAtomCut(1);
  preMicellarAgg_ = 5;
  eulerFlag_ = 0;
  equiMolar_ = 0;
  percolation_ = 0;
}

Space::~Space() {
  checkSizes();
  if (cellType_ > 0) checkCellList();
}

Space* Space::clone() const {
  Space* s = new Space(*this);
  s->reconstruct_();
  return s;
}

shared_ptr<Space> Space::cloneShrPtr() const {
  shared_ptr<Space> s = make_shared<Space>(*this);
  s->reconstruct_();
  return s;
}

void Space::reconstruct_() {
  for (int i = static_cast<int>(addMolList_.size()) - 1; i >= 0; --i) {
    addMolList_[i] = make_shared<Space>(*addMolList_[i]);
  }
  initCellAtomCut(cellAtomCut_);  // this resets ptr neighListChosen_
  Base::reconstruct();
}

int Space::init_config(const int natom) {
  x_.resize(natom*dimen_);
  type_.resize(natom);
  mol_.clear();
  mol2part_.clear();
  nType_.resize(1);
  nType_[0] = natom;
  nMolType_.resize(1);
  nMolType_[0] = natom;

  // low density analytical method of assigning atomic positions
  for (int i = 0; i < natom; ++i) {
    for (int j = 0; j < dimen_; ++j) {
      x_[dimen_*i+j] = 0.95 * (i * dimen_ + j);
    }
    listAtoms_.push_back(i);
    mol_.push_back(i);
    vector<int> m(1, i);
    std::string type("atom");
    moltype_.push_back(type);
    molid_.push_back(0);
    mol2part_.push_back(i);
  }
  mol2part_.push_back(natom);
  for (int i = 0; i < dimen_; ++i) {
    l_[i] = natom * dimen_;
  }

  // update molecule numbers
  xMolGen();

  return 0;
}

int Space::printXYZ(const char* fileName, const int initFlag,
  const std::string comment) {
  stringstream ss;
  ss << fileName << ".xyz";
  FILE * xyzFile = NULL;
  if (initFlag == 1) {
    fileBackUp(ss.str().c_str());
    xyzFile = fopen(ss.str().c_str(), "w");
  } else if (initFlag == 2) {
    if (fileExists(ss.str().c_str())) {
      return 0;
    } else {
      xyzFile = fopen(ss.str().c_str(), "w");
    }
  } else if (initFlag == 0) {
    xyzFile = fopen(ss.str().c_str(), "a");
  } else {
    ASSERT(0, "Unrecognized initFlag:" << initFlag);
  }

  fprintf(xyzFile, "%d\n1 %s\n", natom(), comment.c_str());
  if (xyzFile != NULL) {
    for (int ipart = 0; ipart < natom(); ++ipart) {
      if (type_[ipart] == 0) {
        fprintf(xyzFile, "O ");
      } else if (type_[ipart] == 1) {
        fprintf(xyzFile, "H ");
      } else {
        fprintf(xyzFile, "N ");
      }
      for (int i = 0; i < dimen_; ++i) {
        fprintf(xyzFile, "%f ", x_.at(dimen_*ipart+i));
      }

      // print constant plane for <3D
      for (int i = dimen_; i < 3; ++i) {
        fprintf(xyzFile, "%f ", 0.);
      }
      fprintf(xyzFile, "\n");
    }
  }
  fclose(xyzFile);

  printxyzvmd(fileName, initFlag);
  return 0;
}

void Space::readXYZ(const char* fileName) {
  // check for errors, open if first time
  if (xyzFile_ == NULL) {
    xyzFile_ = make_shared<std::ifstream>(fileName);
    ASSERT(xyzFile_->good(), "cannot open xyz file to read " << fileName);
  }

  // read first two lines
  int iAtom;
  string line;
  getline(*xyzFile_, line);
  if (!xyzFileEOF()) {
    {
      std::istringstream iss(line); iss >> iAtom;
      x_.resize(iAtom*dimen_);
      type_.resize(iAtom);
      nType_.resize(1, iAtom);
      mol_.resize(iAtom);
    }
    getline(*xyzFile_, line);
    { std::istringstream iss(line); iss >> iAtom;}

    // read coordinates and convert to angstroms
    double coord[dimen_];
    listAtoms_.clear();
    int nTypes = 0;
    vector<string> nTypeStrings;
    nType_.clear();
    for (int i = 0; i < natom(); ++i) {
      getline(*xyzFile_, line);
      std::istringstream iss(line);
      string tmp;
      iss >> tmp >> coord[0] >> coord[1] >> coord[2];
      for (int dim = 0; dim < dimen_; ++dim) x_[dimen_*i+dim] = coord[dim];
      listAtoms_.push_back(i);
      if (nTypes == 0) {
        type_[i] = nTypes;
        nType_.push_back(1);
        ++nTypes;
        nTypeStrings.push_back(tmp);
      } else {
        int found = 0;
        for (int t = 0; t < nTypes; ++t) {
          if (tmp.compare(nTypeStrings[t].c_str()) == 0) {
            type_[i] = t;
            ++nType_[t];
            ++found;
          }
        }
        if (found == 0) {
          type_[i] = nTypes;
          nType_.push_back(1);
          ++nTypes;
          nTypeStrings.push_back(tmp);
        } else if (found > 1) {
          ASSERT(0, "more than one(" << found <<
                 ") match in determining type_ in readxyz");
        }
      }
    }

    // update molecule numbers
    xMolGen();
  }
}

void Space::readXYZ(std::ifstream& file) {
  // read first two lines
  int iAtom;
  string line;
  getline(file, line);
  if (!file.eof()) {
    if (natom() == 0) {
      std::istringstream iss(line); iss >> iAtom;

      // add particles
      const double nMol = static_cast<double>(iAtom) /
        static_cast<double>(addMolList_.front()->natom());
      for (int i = 0; i < nMol; ++i) {
        addMol(addMolListType_.front().c_str());
      }
    }
    getline(file, line);
    { std::istringstream iss(line); iss >> iAtom;}

    // read coordinates and convert to angstroms
    double coord[dimen_];
    for (int i = 0; i < natom(); ++i) {
      getline(file, line);
      std::istringstream iss(line);
      string tmp;
      iss >> tmp >> coord[0] >> coord[1] >> coord[2];
      for (int dim = 0; dim < dimen_; ++dim) x_[dimen_*i+dim] = coord[dim];
    }
  }
}

double Space::pbc(const double x,const int i) {
  ASSERT(fabs(xyTilt_)+fabs(xzTilt_)+fabs(yzTilt_) < doubleTolerance,
    "orthogonal box pbc called with nonzero tilt factors");
  double dx = x/l_[i];  // change in position, to be returned
  if (dx > 0.5) {
    dx = -l_[i] * static_cast<int>(dx + 0.5);
  } else if (dx < -0.5) {
    dx = -l_[i] * static_cast<int>(dx - 0.5);
  } else {
    dx = 0.;
  }
  return dx;
}

vector<double> Space::pbc(const vector<double> x) {
  vector<double> dx(dimen_, 0.);
  if (fabs(xyTilt_)+fabs(xzTilt_)+fabs(yzTilt_) < doubleTolerance) {
    for (int dim = 0; dim < dimen_; ++dim) {
      dx[dim] = pbc(x[dim], dim);
    }
  } else {
    vector<double> xnew = x;
    if (dimen_ >= 3) {
      if (fabs(xnew[2]) > 0.5*l_[2]) {
        if (xnew[2] < 0.) {
          dx[2] += l_[2];
          dx[1] += yzTilt_;
          dx[0] += xzTilt_;
          xnew[2] += l_[2];
          xnew[1] += yzTilt_;
          xnew[0] += xzTilt_;
        } else {
          dx[2] -= l_[2];
          dx[1] -= yzTilt_;
          dx[0] -= xzTilt_;
          xnew[2] -= l_[2];
          xnew[1] -= yzTilt_;
          xnew[0] -= xzTilt_;
        }
      }
    }
    if (dimen_ >= 2) {
      if (fabs(xnew[1]) > 0.5*l_[1]) {
        if (xnew[1] < 0.) {
          dx[1] += l_[1];
          dx[0] += xyTilt_;
          xnew[1] += l_[1];
          xnew[0] += xyTilt_;
        } else {
          dx[1] -= l_[1];
          dx[0] -= xyTilt_;
          xnew[1] -= l_[1];
          xnew[0] -= xyTilt_;
        }
      }
    }
    if (dimen_ >= 1) {
      if (fabs(xnew[0]) > 0.5*l_[0]) {
        if (xnew[0] < 0.) {
          dx[0] += l_[0];
          xnew[0] += l_[0];
        } else {
          dx[0] -= l_[0];
          xnew[0] -= l_[0];
        }
      }
    }
  }
//  cout << "dx " << dx[0] << " " << dx[1] << endl;
  return dx;
}

void Space::randDisp(const vector<int> mpart, const double maxDisp) {
  double maxDispTmp = maxDisp;
  for (int dim = 0; dim < dimen_; ++dim) {
    if (maxDisp == -1) maxDispTmp = l_[dim]/2.;
    const double disp = maxDispTmp*(2*uniformRanNum() - 1);
    for (unsigned int i = 0; i < mpart.size(); ++i) {
      x_[dimen_*mpart[i]+dim] += disp;
    }
  }

  // wrap inside box
  wrap(mpart);
}

void Space::randDisp(const int part, const double maxDisp) {
  vector<int> mpart;
  mpart.push_back(part);
  randDisp(mpart, maxDisp);
}

void Space::randDispMulti(const vector<int> mpart, const double maxDisp) {
  double maxDispTmp = maxDisp;
  for (int dim = 0; dim < dimen_; ++dim) {
    if (maxDisp == -1) maxDispTmp = l_[dim]/2.;
    const double disp = maxDispTmp*(2*uniformRanNum() - 1);
    for (unsigned int i = 0; i < mpart.size(); ++i) {
      x_[dimen_*mpart[i]+dim] += disp;
    }
  }
}

void Space::randRotate(const vector<int> mpart, const double maxDisp) {
  if (sphereSymMol_ == false) {
    // assume that mpart is made of only one molecule
    const int iMol = mol_[mpart[0]];

    if (dimen_ == 3) {
      if (eulerFlag_ == 0) {
        // perturb quaternions by factor maxDisp, or,
        // completely randomly, if maxDisp <= 0
        vector<double> q;
        vector<double> qran = quatRandom();
        if (maxDisp > 0) {
          for (int i = 0; i < qdim_; ++i) {
            q.push_back(qMol_[iMol*qdim_+i] + maxDisp*qran[i]);
          }
          const double qsize = sqrt(vecDotProd(q, q));
          for (int i = 0; i < qdim_; ++i) {
            q[i] /= qsize;
          }
        } else {
          q = qran;
        }
        for (int i = 0; i < qdim_; ++i) {
          qMol_[iMol*qdim_+i] = q[i];
        }
      } else {
        // perturb euler angles by a factor maxDisp, or completely randomly
        vector<double> eran = eulerRandom(), e = eran;
        ASSERT(maxDisp <= 0, "euler angle perturbation not implemented. "
          << "Use quaternions instead");
        for (int i = 0; i < qdim_-1; ++i) {
          qMol_[iMol*qdim_+i] = e[i];
        }
      }
    } else if (dimen_ == 2) {
      double theta = qMol_[iMol] + maxDisp*(uniformRanNum() - 0.5);
      theta += pbc2d(theta);
      qMol_[iMol] = theta;
    }

    // update positions with new quaternions
    quat2pos(iMol);
  }
}

void Space::randRotateMulti(const vector<int> mpart,const double maxDisp,
  const vector<double> &sig) {
  const int natom = static_cast<int>(mpart.size());

  // find center of (mass of == 1) mpart
  vector<double> rcm(dimen_, 0.);
  int natomWithMass = 0;
  for (int i = 0; i < natom; ++i) {
    const int iatom = mpart[i];
    for (int dim = 0; dim < dimen_; ++dim) {
      if (static_cast<int>(sig.size()) == 0) {
        rcm[dim] += xcluster_[dimen_*iatom+dim] / static_cast<double>(natom);
      } else {
        if (fabs(sig[type_[iatom]]) > doubleTolerance) {
          ++natomWithMass;
          rcm[dim] += xcluster_[dimen_*iatom+dim];
        }
      }
    }
  }
  for (int dim = 0; dim < dimen_; ++dim) {
    rcm[dim] /= static_cast<double>(natomWithMass)/3.;
  }

  // define all positions as relative to com
  vector<vector<double> > r(natom, vector<double>(dimen_, 0.));
  for (int i = 0; i < natom; ++i) {
    for (int dim = 0; dim < dimen_; ++dim) {
      r[i][dim] = xcluster_[dimen_*mpart[i]+dim] - rcm[dim];
    }
  }

  // find rotation matrix as perturbation
  vector<vector<double> > rot;
  if (dimen_ == 3) {
    vector<double> qran = quatRandom();
    vector<double> qunit(qdim_);
    qunit[qdim_-1] = 1.;
    for (int i = 0; i < qdim_; ++i) {
      qran[i] = qunit[i] + maxDisp*qran[i];
    }
    double qsize = sqrt(vecDotProd(qran, qran));
    for (int i = 0; i < qdim_; ++i) {
      qran[i] /= qsize;
    }
    rot = quat2rot(qran);
  } else if (dimen_ == 2) {
    double theta = maxDisp*(uniformRanNum() - 0.5);
    rot = theta2rot(theta);
  }

  // rotate all positions of mpart about the center of mass
  vector<vector<double> > rnew = matMul(r, rot);
  for (int i = 0; i < natom; ++i) {
    for (int dim = 0; dim < dimen_; ++dim) {
      x_[dimen_*mpart[i]+dim] = rcm[dim] + rnew[i][dim];
    }
  }

  // update molecular orientation
  if (!sphereSymMol_) {
    // are molecules spherically symmetric?
    // no rotations or quaternions necessary
    // set quaternions to unit (0, 0, 0, 1) and xref to x-rcm
    vector<int> molList = mpart2mmol(mpart);

    // obtain the vector from the center of iMol to the point of rotation
    vector<vector<double> > rMol;
    if (eulerFlag_ == 1) {
      int iMolPrev = -1;
      for (unsigned int ipart = 0; ipart < mpart.size(); ++ipart) {
        const int iMol = mol_[mpart[ipart]];
        if (iMol != iMolPrev) {
          rMol.push_back(r[ipart]);
        }
        iMolPrev = iMol;
      }
      ASSERT(rMol.size() == molList.size(), "rMols(" << rMol.size()
             << ")!=molLists(" << molList.size() << ")");
    }

    for (unsigned int i = 0; i < molList.size(); ++i) {
      const int iMol = molList[i];

      // if not using euler angles, define new reference state
      if (eulerFlag_ == 0) {
        if (dimen_ == 3) {
          for (int iq = 0; iq < qdim_ - 1; ++iq) qMol_[qdim_*iMol+iq] = 0;
          qMol_[qdim_*iMol+qdim_-1] = 1;
        } else if (dimen_ == 2) {
          qMol_[iMol] = 0.;
        }
        vector<double> rcmmol(dimen_);
        const vector<int> mmpart = imol2mpart(iMol);
        for (int dim = 0; dim < dimen_; ++dim) {
          rcmmol[dim] = x_[dimen_*mmpart[0]+dim];
        }
        int index = 0;
        for (int ipart = mol2part_[iMol]; ipart < mol2part_[iMol+1]; ++ipart) {
          for (int dim = 0; dim < dimen_; ++dim) {
            xMolRef_[iMol][index][dim] = x_[dimen_*ipart+dim] - rcmmol[dim];
          }
          ++index;
        }
      } else {
        // update euler angles according to the rotation matrix
        vector<vector<double> > Ri = Euler2RotMat(qMol(iMol));
        // invert rotation matrix, because of matMul(r,rot),
        // due to structure of r
        vector<vector<double> > rotT = transpose(rot);
        vector<vector<double> > RiNew = matMul(rotT, Ri);
        vector<vector<double> > euler = RotMat2Euler(RiNew);
        for (int dim = 0; dim < dimen_; ++dim) {
          qMol_[qdim_*iMol+dim] = euler[0][dim];
        }
      }
      quat2pos(iMol);
    }
  }
}

void Space::randRotateMulti(const vector<int> mpart, const double maxDisp) {
  vector<double> sig;
  randRotateMulti(mpart, maxDisp, sig);
}

void Space::delPart(const int ipart) {
  // error check that particle exists
  ASSERT(ipart < natom(), "cannot delete particle that does not exist,"
         << "ipart: " << ipart << " when there are only natom: " << natom());

  // update atom-based cell list
  if ( (fastDel_ == false) && (cellType_ > 0) && (cellAtomCut_) ) {
    eraseAtomFromCell_(ipart);
    atom2cell_.erase(atom2cell_.begin() + ipart);
  }

  // check if ipart is first and only particle on a molecule
  const int iMol = mol_[ipart];
  if ( (fastDel_ == false) &&
       (ipart == mol2part_[iMol]) &&
       (ipart+1 == mol2part_[iMol+1]) ) {
    // delete molecular arrays
//    cout << "deling mol arrays for imol " << iMol << endl;
    moltype_.erase(moltype_.begin() + iMol);
    nMolType_[molid_[iMol]]--;
    molid_.erase(molid_.begin() + iMol);
    if (sphereSymMol_ == false) {
      xMolRef_.erase(xMolRef_.begin() + iMol);
      for (int qd = 0; qd < qdim_; ++qd) {
        qMol_.erase(qMol_.begin() + qdim_*iMol);
      }
    }
    listMols_.erase(listMols_.end() - 1);
    if ( (cellType_ > 0) && (cellAtomCut_ == false) ) {
      // cout << "deling from cell, iMol " << iMol << endl;
      eraseMolFromCell_(iMol);
      mol2cell_.erase(mol2cell_.begin() + iMol);
    }

    // cout down iMol tags in cellList
    for (vector<vector<int> >::iterator it = cellList_.begin();
      it != cellList_.end(); ++it) {
      for (vector<int>::iterator it2 = it->begin(); it2 != it->end(); ++it2) {
        if (*it2 > iMol) {
          // cout << "count down imol " << iMol << " it2 " << *it2 << endl;
          --(*it2);
        }
      }
    }
  }


  // cout down ipart tags in mol2part and tag
  if (fastDel_ == false) {
    for (unsigned int i = 0; i < tag_.size(); ++i) {
      if (tag_[i] == ipart) tag_[i] = -1;
      if (tag_[i] > ipart) --tag_[i];
    }

    for (vector<int>::iterator it = mol2part_.begin(); it != mol2part_.end();
         ++it) {
      if (*it == ipart) {
        mol2part_.erase(it);
      }
    }
    for (vector<int>::iterator it = mol2part_.begin(); it != mol2part_.end();
         ++it) {
      if (*it > ipart) {
        --(*it);
      }
    }
  }

  for (int dim = 0; dim < dimen_; ++dim) x_.erase(x_.begin() + dimen_*ipart);
  --nType_[type_[ipart]];
  type_.erase(type_.begin() + ipart);
  mol_.erase(mol_.begin() + ipart);
  listAtoms_.pop_back();
  mol2part_.back() = natom();

  // count down tags in mol_
  if (fastDel_ == false) {
    int prevMol = -1;
    int mod = 0;
    bool delMolFound = true;
    for (vector<int>::iterator it = mol_.begin(); it != mol_.end(); ++it) {
      if (prevMol == *it) {
        // good
      } else if (prevMol + 1 == *it) {
        ++prevMol;
      } else if (prevMol + 2 == *it) {
        if (delMolFound) {
          int iMol = prevMol;
          if (iMol == -1) iMol = 0;
        }
        delMolFound = false;
        ++prevMol;
        mod = 1;
      } else {
        ASSERT(0, "the structure of mol_ is not continuous when previous mol"
               << "is " << prevMol << " and current mol is " << *it);
      }
      if (mod != 0) *it -= mod;
    }
  }
}

bool Space::fastDelApplicable(const vector<int> mpart) const {
  const int iMol = mol_[mpart.front()];
  if (mol_[mpart.back()] == iMol) {
    if (moltype_[iMol].compare(moltype_.back().c_str()) == 0) {
      return true;  // comment this line to disable fastDel_ across all classes
    }
  }
  return false;
}

void Space::delPart(const vector<int> mpart) {
  fastDel_ = false;
  const int iMol = mol_[mpart.front()];
  if (fastDelApplicable(mpart)) {
    fastDel_ = true;
    fastDelMol_ = iMol;
    const int jMol = nMol() - 1;

    // swap mpart with last molecule
    for (int i = static_cast<int>(mpart.size()) - 1; i >= 0; --i) {
      const int ipart = mpart[i];
      const int jpart = natom() - static_cast<int>(mpart.size()) + i;
      for (int dim = 0; dim < dimen_; ++dim) {
        x_[dimen_*ipart+dim] = x_[dimen_*jpart+dim];
      }
      const int t = type_[ipart];
      type_[ipart] = type_[jpart];
      type_[jpart] = t;

      // update tag if deleted
      for (unsigned int i = 0; i < tag_.size(); ++i) {
        if (tag_[i] == ipart) tag_[i] = -1;
      }
    }

    // update cell list
    // cout << "there are " << nMol() << " mols currently " << endl;
    if (cellType_ > 0) {
      // cout << "fastdel updating cell list: for iMol" << iMol << endl;
      updateCellofiMol(iMol);
      if (cellAtomCut_) {
        for (int i = static_cast<int>(mpart.size()) - 1; i >= 0; --i) {
          eraseAtomFromCell_(natom() - static_cast<int>(mpart.size()) + i);
        }
        for (int i = static_cast<int>(mpart.size()) - 1; i >= 0; --i) {
          atom2cell_.pop_back();
        }
      } else {
        // cout << "erasing  nMol()-1= " << nMol() - 1 << endl;
        eraseMolFromCell_(nMol() - 1);
        mol2cell_.pop_back();
      }
    }

    if (sphereSymMol_ == false) {
      // cout << "deling qMol and xMol for iMol "
      // << iMol << " jMol " << jMol << endl;
      // cout << "qMol beginning size " << qMol_.size() << endl;
      for (int qd = 0; qd < qdim_; ++qd) {
        qMol_[qdim_*iMol+qd] = qMol_[qdim_*jMol+qd];
      }
      for (int qd = 0; qd < qdim_; ++qd) {
        qMol_.pop_back();
      }
      xMolRef_[iMol] = xMolRef_[jMol];
      xMolRef_.pop_back();
    }
    mol2part_.pop_back();
    moltype_.pop_back();
    nMolType_[molid_.back()]--;
    molid_.pop_back();
    listMols_.pop_back();

    // delete last molecule without any countdown of tags
    for (int i = static_cast<int>(mpart.size()) - 1; i >= 0; --i) {
      delPart(natom()-1);
    }
  }
  if (fastDel_ == false) {
    for (int i = static_cast<int>(mpart.size()) - 1; i >= 0; --i) {
      delPart(mpart[i]);
    }
  }
}

void Space::addPart(const vector<double> v, const int itype, const int imol) {
  // error check that dimensions match
  ASSERT(static_cast<int>(v.size()) == dimen_, "dimensions of space ("
         << dimen_ << ") and position of new particle ("
         << v.size() << ") do not match");

  for (int dim = 0; dim < dimen_; ++dim) x_.push_back(v[dim]);
  type_.push_back(itype);
  if (itype > nParticleTypes() - 1) nType_.resize(itype + 1);
  ++nType_[itype];
  mol_.push_back(imol);
  listAtoms_.push_back(natom() - 1);
}

void Space::readXYZBulk(const int nMolAtoms, const char* type,
  const char* fileName) {
  readxyz(fileName);
  std::string typestr(type);
  int molid = -1;
  vector<int> m;
  mol2part_.clear();
  nType_.clear();
  nType_.resize(2, 0);
  for (int i = 0; i < natom(); i++) {
    if (i % nMolAtoms == 0) {
      type_.at(i) = 0;
      ++nType_[0];
      ++molid;
      mol2part_.push_back(i);
    } else {
      type_.at(i) = 1;
      ++nType_[1];
    }
    mol_.at(i) = molid;
    if (i % nMolAtoms == nMolAtoms - 1) {
      moltype_.push_back(typestr);
      molid_.push_back(0);    // place holder
    }
  }
  mol2part_.push_back(natom());

  // update molecule numbers, xMol, quaternions and molecule reference vectors
  qMolInit();
}

vector<int> Space::randMol() {
  vector<int> mpart;
  const int iMol = uniformRanNum(0, nMol() - 1);
  for (int i = mol2part_[iMol]; i < mol2part_[iMol+1]; ++i) mpart.push_back(i);
  return mpart;
}

vector<int> Space::randMolDiff(const vector<int> jmMol) {
  vector<int> mpart;

  ASSERT(static_cast<int>(jmMol.size()) != nMol(),
         "cannot pick random molecule while excluding all molecules");

  // pick a random mol number from random atom
  // check that the molecule is not part of the same molecular
  // as any particles in jmpart
  int term = 0;
  while (term == 0) {
    mpart = randMol();
    const int imol = mol_[mpart[0]];
    term = 1;
    for (unsigned int i = 0; i < jmMol.size(); ++i) {
      if (jmMol[i] == imol) {
        term = 0;
      }
    }
  }

  return mpart;
}

vector<int> Space::randMolDiff(const int jMol) {
  const vector<int> jMolVec(1, jMol);
  return randMolDiff(jMolVec);
}

vector<int> Space::randMolSubset(const vector<int> jmMol) {
  ASSERT(static_cast<int>(jmMol.size()) != 0,
         "cannot pick random molecule from null subset of molecules");
  vector<int> mpart;
  const int jrandMol = jmMol[uniformRanNum(0,
                       static_cast<int>(jmMol.size()) - 1)];
  for (int ipart = mol2part_[jrandMol]; ipart < mol2part_[jrandMol+1];
         ++ipart) mpart.push_back(ipart);
  return mpart;
}

void Space::restore(vector<int> mpart) {
  ASSERT(xold_.size() >= mpart.size(), "restore() requires previous use of "
    << "xStore()");
  for (unsigned int i = 0; i < mpart.size(); ++i) {
    for (int dim = 0; dim < dimen_; ++dim) {
      x_[dimen_*mpart[i]+dim] = xold_[i][dim];
    }
  }
  if (sphereSymMol_ == false) {
    // assume that mpart is made of only one molecule
    const int iMol = mol_[mpart[0]];
    for (int dim = 0; dim < qdim_; ++dim) {
      qMol_[qdim_*iMol+dim] = qMolOld_[0][dim];
    }
  }
}

void Space::restoreAll() {
  // NOTE HWH ASSERT for size
  ASSERT(xOldAll_.size() == x_.size(), "stored particle coordinates size "
    << xOldAll_.size() << " does not match current size " << x_.size());
  x_ = xOldAll_;
  if (sphereSymMol_ == false) {
    ASSERT(qMol_.size() == qMolOldAll_.size(), "size mismatch");
    qMol_ = qMolOldAll_;
    ASSERT(xMolRef_.size() == xMolRefOld_.size(), "size mismatch");
    xMolRef_ = xMolRefOld_;
  }
}

void Space::addMol(const char* type) {
  std::string typestr(type);
  vector<vector<double> > xmol;

  shared_ptr<Space> s = findAddMolInList(typestr);

  if (sphereSymMol_ && (s->natom() > 1)) sphereSymMol_ = false;
  vector<vector<double> > xn = s->xMol().front();
  if (!sphereSymMol_) xMolRef_.push_back(s->xMol().front());

  // move molecule to random position within the box
  //  or if xAdd !null, to position xadd
  for (int dim = 0; dim < dimen_; ++dim) {
    double disp;
    if (static_cast<int>(xAdd.size()) == 0) {
      disp = l_[dim]*(uniformRanNum() - 0.5);
    } else {
      disp = xAdd[dim];
    }
    for (unsigned int ipart = 0; ipart < xn.size(); ++ipart) {
      xn[ipart][dim] += disp;
    }
  }

  // add particles
  int imol = 0;
  if (natom() > 0) imol = mol_.back() + 1;
  for (unsigned int ipart = 0; ipart < xn.size(); ++ipart) {
    addPart(xn[ipart], s->type()[ipart], imol);
  }
  moltype_.push_back(s->moltype().front());
  molid_.push_back(findAddMolListIndex(typestr));

  // update number of molecules
  mol2part_.push_back(natom());
  listMols_.push_back(nMol() - 1);
  if (molid_.back() > static_cast<int>(nMolType_.size())-1) {
    nMolType_.resize(molid_.back() + 1);
  }
  ++nMolType_[molid_.back()];

  // reference quaternion
  if (sphereSymMol_ == false) {
    if (static_cast<int>(xAdd.size()) == 0) {
      if (dimen_ == 3) {
        if (eulerFlag_ == 0) {
          vector<double> qran = quatRandom();
          for (int qd = 0; qd < qdim_; ++qd) qMol_.push_back(qran[qd]);
        } else {
          vector<double> eran = eulerRandom();
          for (int qd = 0; qd < qdim_-1; ++qd) qMol_.push_back(eran[qd]);
          qMol_.push_back(1);
        }
      } else if (dimen_ == 2) {
        qMol_.push_back(2*PI*uniformRanNum());
      }
    } else {
      if (dimen_ == 3) {
        for (int qd = 0; qd < qdim_ - 1; ++qd) qMol_.push_back(0);
        qMol_.push_back(1);
      } else if (dimen_ == 2) {
        qMol_.push_back(0);
      }
    }
    quat2pos(nMol()-1);
  }

  // update cell list
  if (cellType_ > 0) {
    if (cellAtomCut_) {
      for (int ipart = mol2part_[nMol()-1]; ipart < mol2part_[nMol()];
           ++ipart) {
        addAtomtoCell_(ipart);
      }
    } else {
      addMoltoCell_(nMol() - 1);
    }
  }

  // xAdd is used to add a particle to a specific location
  // once that is accomplished, clear it for random add
  xAdd.clear();

  // wrap in box
  const vector<int> mpart = lastMolIDVec();
  wrap(mpart);
}

/**
 * initialize the potential addition of molecules
 */
void Space::addMolInit(const char* fileName   //!< data file for molecule
  ) {
  ASSERT(fileExists(fileName),
         "addMolInit(fileName = " << fileName << ") doesn't exists.");
  addMolList_.push_back(std::make_shared<Space>(dimen_, 0));

  // determine the number of particle types that already exist,
  //  thus, creating new particles with types higher than those existing
  int nTypesExist = 0;
  if (addMolList_.size() > 1) nTypesExist = nParticleTypes();

  addMolList_.back()->initData(fileName, nTypesExist);
  std::string fs(fileName);
  addMolListType_.push_back(fs);

  // add and delete molecule to initialize nType arrays
  addMol(fileName);
  delPart(lastMolIDVec());
}

void Space::xMolGen() {
  xMol_.clear();
  if (natom() != 0) {
    vector<vector<double> > xMolOne;
    int molPrev = -1;
    molPrev = mol_.front();
    for (int ipart = 0; ipart < natom(); ++ipart) {
      if (molPrev != mol_[ipart]) {
        xMol_.push_back(xMolOne);
        xMolOne.clear();
      }
      vector<double> xVec;
      for (int dim = 0; dim < dimen_; ++dim) {
        xVec.push_back(x_[dimen_*ipart+dim]);
      }
      xMolOne.push_back(xVec);
      molPrev = mol_[ipart];
      if (ipart == natom() -1) xMol_.push_back(xMolOne);
    }
  }
  listMols_.clear();
  for (int i = 0; i < nMol(); ++i) listMols_.push_back(i);
}

int Space::checkBond(const char* type, const double tol) {
  std::string typestr(type);

  vector<vector<double> > xRef;
  if (typestr.compare("spce") == 0) {
    xRef = vecSPCE();
  } else {
    ASSERT(0, "unrecognized molecule type for checkBond function of space.cc");
  }

  // start as bonds match, and switch to zero if a test fails
  int bondsMatch = 1;

  // generate molecule list
  xMolGen();

  // searching all molecules of same size as xRef, check bond lengths
  for (int imol = 0; imol < nMol(); ++imol) {
    const int nPart = static_cast<int>(xRef.size());
    if (static_cast<int>(xMol_[imol].size()) == nPart) {
      for (int i = 0; i < nPart; ++i) {
        for (int j = 0; j < nPart; ++j) {
          double r, rref, r2 = 0., r2ref = 0.;
          for (int dim = 0; dim < dimen_; ++dim) {
            r = xMol_[imol][i][dim] - xMol_[imol][j][dim];
            rref = xRef[i][dim] - xRef[j][dim];
            r2 += r*r;
            r2ref += rref*rref;
          }
          if (fabs(r2 - r2ref) > tol) {
            cout << "r2 r2ref " << r2 << " " << r2ref << " fabs "
                 << fabs(r2 - r2ref) << endl;
            ASSERT(0, "checkBond failed");
            bondsMatch = 0;
          }
        }
      }
    }
  }

  // check that all quaternions are normalized
  if ( (eulerFlag_ == 0) && (static_cast<int>(qMol_.size()) != 0) ) {
    for (int iMol = 0; iMol < nMol(); ++iMol) {
      double q2 = 0;
      for (int iq = 0; iq < qdim_; ++iq) {
        q2 += qMol(iMol, iq);
      }
      if (fabs(q2-1) > tol) {
        cout << "q2(" << q2 << ") != 1" << endl;
        ASSERT(0, "checkBond failed");
        bondsMatch = 0;
      }
    }
  }

  return bondsMatch;
}

int Space::checkBond(const double tol) {
  // start as bonds match, and switch to zero if a test fails
  int bondsMatch = 1;
  if (nMol() > 0) {
    // generate molecule list
    xMolGen();

    // searching all molecules of same size as xRef, check bond lengths
    for (int imol = 0; imol < nMol(); ++imol) {
      shared_ptr<Space> s = findAddMolInList(moltype_[imol]);
      for (unsigned int i = 0; i < xMol_[imol].size(); ++i) {
        for (int j = 0; j < s->natom(); ++j) {
          double r, rref, r2 = 0., r2ref = 0.;
          for (int dim = 0; dim < dimen_; ++dim) {
            r = xMol_[imol][i][dim] - xMol_[imol][j][dim];
            rref = s->xMol()[0][i][dim] - s->xMol()[0][j][dim];
            r2 += r*r;
            r2ref += rref*rref;
          }
          if (fabs(r2 - r2ref) > tol) {
            ASSERT(0, "r2 r2ref " << r2 << " " << r2ref << " fabs "
              << fabs(r2 - r2ref)
              << "checkBond failed between atom " << i << " and atom " << j);
            bondsMatch = 0;
          }
        }
      }
    }
  }
  return bondsMatch;
}

double Space::minl() const {
  double minl = 1e15;
  for (int dim = 0; dim < dimen_; ++dim) {
    minl = std::min(minl, l_[dim]);
  }
  return minl;
}

int Space::bonded(const int iAtom, const int jAtom,
                  const double rAbove, const double rBelow) {
  vector<double> xm, xj;
  for (int dim = 0; dim < dimen_; ++dim) {
    xm.push_back(x_[dimen_*iAtom+dim]);
    xj.push_back(x_[dimen_*jAtom+dim]);
  }
  const double r2 = rsq(xm, xj);
  if ( (r2 < rAbove*rAbove) && (r2 > rBelow*rBelow) ) {
    return 1;
  } else {
    return 0;
  }
}

int Space::bonded(const vector<int> mpart, const vector<int> jmpart,
                  const double rabove, const double rbelow) {
  return bonded(mpart[0], jmpart[0], rabove, rbelow);
}

int Space::avb(const vector<int> mpart, const vector<int> jmpart,
               const double rabove, const double rbelow, const char* type) {
  ASSERT(mpart[0] != jmpart[0], "Aggregation volume bias move attempted"
    << "where the first atom in the moving particle, mpart(" << mpart[0]
    << ") is equilvalent to the first atom in the target particle, jmpart("
    << jmpart[0] << ").");

  // find whether mpart and jpart are bonded before move
  int bond = bonded(mpart, jmpart, rabove, rbelow);

  // store old position
  vector<double> xold;
  for (int dim = 0; dim < dimen_; ++dim) {
    xold.push_back(x_[dimen_*mpart[0]+dim]);
  }

  // compute new position (xnew) of first particle in molecule mpart
  // to move to bonded or unbonded region of jmpart
  std::string typestr(type);
  vector<double> xnew;
  if (typestr.compare("bonded") == 0) {
    // move mpart to bonded position
    xnew = ranShell(rabove, rbelow, dimen_);
    for (int dim = 0; dim < dimen_; ++dim) {
      xnew[dim] += x_[dimen_*jmpart[0]+dim];
    }

  } else if (typestr.compare("nonbonded") == 0) {
    // move mpart to nonbonded position
    // try any position in box, reject if it ends up in bonded region
    int term = 0;
    xnew.resize(dimen_, 0.);
    vector<double> xj;
    for (int dim = 0; dim < dimen_; ++dim) {
      xj.push_back(x_[dimen_*jmpart[0]+dim]);
    }
    while (term == 0) {
      for (int dim = 0; dim < dimen_; ++dim) {
        xnew[dim] = l_[dim]*(uniformRanNum() - 0.5);
      }
      const double r2 = rsq(xnew, xj);
      if ( (r2 > rabove*rabove) || (r2 < rbelow*rbelow) ) {
        term = 1;
      }
    }
  } else {
    ASSERT(0, "unrecognized type " << type << " in avb space.cc");
  }

  // rigidly translate mpart such that first particle is positioned at xnew
  for (unsigned int i = 0; i < mpart.size(); ++i) {
    const int ipart = mpart[i];
    for (int dim = 0; dim < dimen_; ++dim) {
      x_[dimen_*ipart+dim] += xnew[dim] - xold[dim];
    }
  }

  // randomly orient mpart
  randRotate(mpart, -1);

  // wrap new particle position inside box
  wrap(mpart);

  return bond;
}

void Space::xStore(const vector<int> mpart   //!< list of particles to store
  ) {
  xold_.resize(static_cast<int>(mpart.size()), vector<double>(dimen()));
  for (unsigned int i = 0; i < mpart.size(); ++i) {
    for (int dim = 0; dim < dimen_; ++dim) {
      xold_[i][dim] = x_[dimen_*mpart[i]+dim];
    }
  }
  if (sphereSymMol_ == false) {
    // assume that mpart is made of only one molecule
    const int iMol = mol_[mpart[0]];
    qMolOld_.resize(1, vector<double>(qdim()));
    for (int dim = 0; dim < qdim_; ++dim) {
      qMolOld_[0][dim] = qMol_[qdim_*iMol+dim];
    }
  }
}

void Space::xStoreAll() {
  xOldAll_ = x_;
  if (sphereSymMol_ == false) {
    qMolOldAll_ = qMol_;
    xMolRefOld_ = xMolRef_;
  }
}

void Space::xStoreMulti(const vector<int> mpart, const int flag) {
  // restore if flag is positive
  if (flag >= 0) {
    ASSERT(mpart.size() == xOldMulti_[flag].size(), "size mismatch");
    for (unsigned int i = 0; i < mpart.size(); ++i) {
      for (int dim = 0; dim < dimen_; ++dim) {
        x_[dimen_*mpart[i]+dim] = xOldMulti_[flag][i][dim];
      }
    }
    if (sphereSymMol_ == false) {
      // assume that mpart is made of only one molecule
      const int iMol = mol_[mpart[0]];
      for (int dim = 0; dim < qdim_; ++dim) {
        qMol_[qdim_*iMol+dim] = qMolOldMulti_[flag][0][dim];
      }
    }

  // store if flag is negative
  } else {
    if (flag == -1) xOldMulti_.clear();
    const int m = static_cast<int>(xOldMulti_.size());
    xOldMulti_.resize(m+1, vector<vector<double> >
      (static_cast<int>(mpart.size()), vector<double>(dimen())));
    for (unsigned int i = 0; i < mpart.size(); ++i) {
      for (int dim = 0; dim < dimen_; ++dim) {
        xOldMulti_[m][i][dim] = x_[dimen_*mpart[i]+dim];
      }
    }
    if (sphereSymMol_ == false) {
      // assume that mpart is made of only one molecule
      const int iMol = mol_[mpart[0]];
      if (flag == -1) qMolOldMulti_.clear();
      const int m = static_cast<int>(qMolOldMulti_.size());
      qMolOldMulti_.resize(m+1, vector<vector<double> >
        (1, vector<double>(qdim())));
      for (int dim = 0; dim < qdim_; ++dim) {
        qMolOldMulti_[m][0][dim] = qMol_[qdim_*iMol+dim];
      }
    }
  }
}

double Space::rsq(const vector<double> xi, const vector<double> xj) {
  vector<double> r(dimen_);
  for (int dim = 0; dim < dimen_; ++dim) {
    r[dim] = xi[dim] - xj[dim];
  }
  const vector<double> dx = pbc(r);
  for (int dim = 0; dim < dimen_; ++dim) {
    r[dim] += dx[dim];
  }
  return vecDotProd(r, r);
}

void Space::wrap(const vector<int> mpart) {
  // find displacement of first particle in molecule in order to wrap
  vector<double> x;
  for (int dim = 0; dim < dimen_; ++dim) {
    x.push_back(x_[dimen_*mpart[0]+dim]);
  }
  vector<double> xold = x;
  rwrap(&x);

  // wrap all particles by this displacement
  for (unsigned int i = 0; i < mpart.size(); ++i) {
    const int ipart = mpart[i];
    for (int dim = 0; dim < dimen_; ++dim) {
      x_[dimen_*ipart+dim] += x[dim] - xold[dim];
    }
  }
  // if (cellType_ > 0) updateCellofiMol(mol_[mpart.front()]);
}

void Space::rwrap(vector<double> *rvecPtr) {
  vector<double>& rvec = *rvecPtr;
  vector<double> dx(dimen_, 0.);
  if (fabs(xyTilt_)+fabs(xzTilt_)+fabs(yzTilt_) < doubleTolerance) {
    dx = pbc(rvec);
    for (int dim = 0; dim < dimen_; ++dim) {
      rvec[dim] += dx[dim];
    }
  } else {
    if (dimen_ == 2) {
      if (rvec[0] > (0.5*l_[0] + rvec[1]/l_[1]*xyTilt_)) {
        dx[0] -= l_[0];
      }
      if (rvec[0] < (-0.5*l_[0] + rvec[1]/l_[1]*xyTilt_)) {
        dx[0] += l_[0];
      }
      if (rvec[1] >  0.5*l_[1]) dx[1] -= l_[1];
      if (rvec[1] < -0.5*l_[1]) dx[1] += l_[1];
      for (int dim = 0; dim < dimen_; ++dim) {
        rvec[dim] += dx[dim];
      }

    } else if (dimen_ == 3) {
      if (rvec[2] >  0.5*l_[2]) {
        dx[2] -= l_[2];
        rvec[2] -= l_[2];
      } else if (rvec[2] < -0.5*l_[2]) {
        dx[2] += l_[2];
        rvec[2] += l_[2];
      }
      if (rvec[1] >  0.5*l_[1] + rvec[2]/l_[2]*yzTilt_) {
        dx[1] -= l_[1];
        rvec[1] -= l_[1];
      } else if (rvec[1] < -0.5*l_[1] + rvec[2]/l_[2]*yzTilt_) {
        dx[1] += l_[1];
        rvec[1] += l_[1];
      }
      if (rvec[0] > (0.5*l_[0] + rvec[1]/l_[1]*xyTilt_
          + rvec[2]/l_[2]*xzTilt_)) {
        dx[0] -= l_[0];
        rvec[0] -= l_[0];
      } else if (rvec[0] < (-0.5*l_[0] + rvec[1]/l_[1]*xyTilt_
                 + rvec[2]/l_[2]*xzTilt_)) {
        dx[0] += l_[0];
        rvec[0] += l_[0];
      }
    }
  }
  // if (fabs(dx[0]) + fabs(dx[1]) + fabs(dx[2]) > 0) {
  //   cout << "WRAPPED~" << endl;
  // }
}

void Space::updateCells(const double dCellMin, const double rCut) {
  if (dCellMin >= rCut) {
    cellType_ = 1;
  } else {
    cellType_ = 2;
  }
  dCellMin_ = dCellMin;
  nCellVec_.clear();
  dCell_.clear();
  for (int dim = 0; dim < dimen_; ++dim) {
    nCellVec_.push_back(
      static_cast<int>(l_[dim]/static_cast<double>(dCellMin)));
    dCell_.push_back(l_[dim] / static_cast<double>(nCellVec_[dim]));
    if ( (dCell_[dim] + rCut > l_[dim]/2.) || (l_[dim] == 0) ) {
      cellType_ = 0;
      WARN(verbose_ == 1, "cell list disabled due to dCell (" << dCell_[dim]
        << ") + rCut (" << rCut << ") > l_[dim]/2 (l_[" << dim << "]="
        << l_[dim] << ")");
    }
  }
  nCell_ = product(nCellVec_);
  if ( (dimen_ == 3) && (cellType_ != 0) ) {
    if (cellType_ == 1) {
      int mix, miy, miz, mjx, mjy, mjz;
      neighCell_.clear();
      neighCell_.resize(nCell_);
      for (mix = 0; mix < nCellVec_[0]; ++mix) {
      for (miy = 0; miy < nCellVec_[1]; ++miy) {
      for (miz = 0; miz < nCellVec_[2]; ++miz) {
        const int icell = mvec2m3d_(mix, miy, miz);
        for (mjx = mix-1; mjx <= mix+1; ++mjx) {
        for (mjy = miy-1; mjy <= miy+1; ++mjy) {
        for (mjz = miz-1; mjz <= miz+1; ++mjz) {
          neighCell_[icell].push_back(mvec2m3d_(mjx, mjy, mjz));
        }}}
      }}}
    } else if (cellType_ == 2) {
      // cout << "enumerate box corners of all cells" << endl;
      vector<vector<vector<double> > > cellCorners;
      for (int m = 0; m < nCell_; ++m) {
        vector<vector<double> > r = cellCorners_(m);
        cellCorners.push_back(r);
      }

      // cout << "looping through all pairs of cells, put neighbors"
      //      << "in neighCell by checking box corners" << endl;
      vector<vector<int> > neighCell(nCell_, vector<int>(nCell_));
      const int nCorners = static_cast<int>(pow(2, dimen_));
      for (int mi = 0; mi < nCell_ - 1; ++mi) {
        for (int mj = mi + 1; mj < nCell_; ++mj) {
          int ci = 0, neigh = 0;
          while ( (neigh == 0) && (ci < nCorners) ) {
            int cj = 0;
            while ( (neigh == 0) && (cj < nCorners) ) {
              double r, r2 = 0.;
              for (int dim = 0; dim < dimen_; ++dim) {
                r = cellCorners[mi][ci][dim] - cellCorners[mj][cj][dim];
                if (r >  0.5 * l_[dim]) r -= l_[dim];
                if (r < -0.5 * l_[dim]) r += l_[dim];
                r2 += r*r;
              }
              if (r2 < rCut*rCut) neigh = 1;
              ++cj;
            }
            ++ci;
          }
          if (neigh == 1) {
            neighCell[mi][mj] = 1;
            neighCell[mj][mi] = 1;
          }
        }
      }
      for (int mi = 0; mi < nCell_; ++mi) neighCell[mi][mi] = 1;

      // build the cell mask pointer array
      cMaskPnt_.clear();
      cMaskPnt_.resize(nCell_);
      for (int mi = 0; mi < 1; ++mi) {
        bool stripeOn = false;
        vector<int> cMaskTmp;
        int nneigh = 0;
        for (int mj = 0; mj < nCell_; ++mj) {
          if (neighCell[mi][mj] == 1) {
            if (stripeOn == false) {
              cMaskTmp.push_back(mj);
              stripeOn = true;
            }
            ++nneigh;
          } else {
            if (stripeOn == true) {
              cMaskTmp.push_back(mj);
              stripeOn = false;
            }
          }
        }
        cMaskPnt_.push_back(cMaskTmp);
      }
    } else {
      ASSERT(0, "Unrecognized cellType_");
    }

    buildCellList();
  } else if ( (dimen_ == 2) && (cellType_ == 1) ) {
    int mix, miy, mjx, mjy;
    neighCell_.clear();
    neighCell_.resize(nCell_);
    for (mix = 0; mix < nCellVec_[0]; ++mix) {
    for (miy = 0; miy < nCellVec_[1]; ++miy) {
      const int icell = mvec2m2d_(mix, miy);
      for (mjx = mix-1; mjx <= mix+1; ++mjx) {
      for (mjy = miy-1; mjy <= miy+1; ++mjy) {
        neighCell_[icell].push_back(mvec2m2d_(mjx, mjy));
      }}
    }}
    buildCellList();
  } else if (cellType_ > 0) {
    ASSERT(0, "Cell list not implemented for this dimensionality");
    cellType_ = 0;
  }
}

vector<int> Space::m2vec_(const int m) {
  vector<int> mVec;
  mVec.push_back(m % nCellVec_[0]);
  if (dimen_ > 1) {
    mVec.push_back(static_cast<int>((m % (nCellVec_[0]*nCellVec_[1]))/
                   nCellVec_[0]));
  }
  if (dimen_ > 2) {
    mVec.push_back(static_cast<int>(m/(nCellVec_[0]*nCellVec_[1])));
  }
  return mVec;
}

vector<vector<double> > Space::cellCorners_(const int m) {
  vector<vector<double> > rm;
  vector<int> mVec = m2vec_(m);
  vector<double> rc;
  for (int dim = 0; dim < dimen_; ++dim) {
    rc.push_back(static_cast<double>
      (mVec[dim] + 0.5)*dCell_[dim] - l_[dim]/2.);
  }
  double cx, cy, cz;
  for (int xi = -1; xi <= 1; xi += 2) {
    cx = static_cast<double>(xi)/2*dCell_[0];
  for (int yi = -1; yi <= 1; yi += 2) {
    cy = static_cast<double>(yi)/2*dCell_[1];
  for (int zi = -1; zi <= 1; zi += 2) {
    cz = static_cast<double>(zi)/2*dCell_[2];
    vector<double> r;
    r.push_back(rc[0] + cx);
    r.push_back(rc[1] + cy);
    r.push_back(rc[2] + cz);
    rm.push_back(r);
  }}}
  return rm;
}

int Space::rvec2m_(const vector<double> &r) {
  int cell = -1;
  if (dimen_ == 3) {
    const int xc = static_cast<int>
      (nCellVec_[0]*(r[0]/l_[0] + 0.5)) % nCellVec_[0];
    const int yc = static_cast<int>
      (nCellVec_[1]*(r[1]/l_[1] + 0.5)) % nCellVec_[1];
    const int zc = static_cast<int>
      (nCellVec_[2]*(r[2]/l_[2] + 0.5)) % nCellVec_[2];
    cell = nCellVec_[0]*nCellVec_[1]*zc +
                     nCellVec_[0]*yc + xc;
    if ( (cell < 0) || (cell > nCell_ - 1) ) {
      ASSERT(0, "cell(" << cell << ") is either <0 or > nCell(" << nCell_
        << " - 1) " << endl << " r " << r[0] << " " << r[1] << " " << r[2]
        << endl << " xc " << xc << " " << yc << " " << zc << endl <<
        " l " << l_[0] << " " << l_[1] << " " << l_[2] << endl << " nCellVec "
        << nCellVec_[0] << " " << nCellVec_[1] << " " << nCellVec_[2]);
    }
  } else if (dimen_ == 2) {
    const int xc = static_cast<int>
      (nCellVec_[0]*(r[0]/l_[0] + 0.5)) % nCellVec_[0];
    const int yc = static_cast<int>
      (nCellVec_[1]*(r[1]/l_[1] + 0.5)) % nCellVec_[1];
    cell = nCellVec_[0]*yc + xc;
    if ( (cell < 0) || (cell > nCell_ - 1) ) {
      ASSERT(0, "cell(" << cell << ") is either <0 or > nCell(" << nCell_
      << " - 1) " << endl <<
        " r " << r[0] << " " << r[1] << endl <<
        " xc " << xc << " " << yc << endl <<
        " l " << l_[0] << " " << l_[1] << endl <<
        " nCellVec " << nCellVec_[0] << " " << nCellVec_[1]);
    }
  }
  return cell;
}

int Space::iatom2m(const double &ipart) {
  vector<double> rvec(dimen_);
  for (int dim = 0; dim < dimen_; ++dim) rvec[dim] = x(ipart, dim);
  rwrap(&rvec);
  return rvec2m_(rvec);
}

void Space::buildCellList() {
  if (cellType_ != 1) {
    ASSERT(0, "cellType other than 1 isn't implemented");
  } else if (cellType_ == 1) {
    cellList_.clear();
    cellList_.resize(nCell_);
    xMolGen();
    if (cellAtomCut_) {
      atom2cell_.clear();
      for (int ipart = 0; ipart < natom(); ++ipart) {
        addAtomtoCell_(ipart);
      }
    } else {
      mol2cell_.clear();
      for (int iMol = 0; iMol < nMol(); ++iMol) {
        addMoltoCell_(iMol);
        // cout << "iMol " << iMol << " x " << xMol_[iMol][0][0] << " "
        //      << xMol_[iMol][0][1] << " " << xMol_[iMol][0][2]
        //      << " iCell " << iCell << endl;
      }
    }
  }
}

vector<int> Space::lastMolIDVec() {
  vector<int> mpart;
  for (int ipart = mol2part_[nMol() - 1]; ipart < mol2part_.back(); ++ipart) {
    mpart.push_back(ipart);
  }
  return mpart;
}

void Space::qMolInit() {
  // update molecule numbers and xMol
  xMolGen();
  xMolRef_ = xMol_;
  sphereSymMol_ = false;
  qMol_.resize(nMol()*qdim_);
  for (int i = 0; i < nMol(); ++i) {
    qMolInit(i);
  }
}

void Space::qMolInit(const int iMol) {
  if (dimen_ == 3) {
    for (int dim = 0; dim < dimen_; ++dim) qMol_[qdim_*iMol+dim] = 0;
    qMol_[qdim_*iMol+dimen_] = 1;
  } else if (dimen_ == 2) {
    qMol_[iMol] = 0;
  }

  // initialize xMolRef_ with positions
  for (int ipart = mol2part_[iMol]; ipart < mol2part_[iMol+1]; ++ipart) {
    const int i = ipart - mol2part_[iMol];
    for (int dim = 0; dim < dimen_; ++dim) {
      xMolRef_[iMol][i][dim] = x(ipart, dim);
    }
  }

  // shift xMolRef_ such that the first atom is zero, e.g. xMolRef[][0][] == 0
  if (xMolRef_[iMol].size() > 1) {
    for (unsigned int i = 1; i < xMolRef_[iMol].size(); ++i) {
      for (int dim = 0; dim < dimen_; ++dim) {
        xMolRef_[iMol][i][dim] -= xMolRef_[iMol][0][dim];
      }
    }
  }
  for (int dim = 0; dim < dimen_; ++dim) xMolRef_[iMol][0][dim] = 0.;
}

void Space::quat2pos(const int iMol) {
  vector<vector<double> > xref = xMolRef_[iMol];
  if (static_cast<int>(xref.size()) > 1) {
    vector<vector<double> > r, xnew;
    if (dimen_ == 3) {
      if (eulerFlag_ == 0) {
        vector<double> q(qdim_);
        for (int i = 0; i < qdim_; ++i) q[i] = qMol_[iMol*qdim_+i];
        r = quat2rot(q);
        xnew = matMul(xref, r);
      } else {
        vector<double> e(dimen_);
        for (int i = 0; i < dimen_; ++i) e[i] = qMol_[iMol*qdim_+i];
        r = Euler2RotMat(e);
        vector<vector<double> > xreftp = transpose(xref);
        xnew = matMul(r, xreftp);
        xnew = transpose(xnew);
      }

    } else if (dimen_ == 2) {
      r = theta2rot(qMol_[iMol]);
      xnew = matMul(xref, r);
    }
    const int iPartPivot = mol2part_[iMol];
    for (unsigned int i = 1; i < xref.size(); ++i) {
      const int ipart = mol2part_[iMol] + i;
      for (int dim = 0; dim < dimen_; ++dim) {
        x_[dimen_*ipart+dim] = x_[dimen_*iPartPivot+dim] + xnew[i][dim];
      }
    }
  }
}

void Space::settype(const int iatom, const int itype) {
  --nType_.at(type_[iatom]);
  if (nParticleTypes() -1 < itype) nType_.resize(itype+1);
  type_.at(iatom) = itype;
  ++nType_.at(itype);
}

void Space::buildNeighListCellAtomCut(const int ipart) {
  neighListCell_.clear();
  neighListChosen_ = &neighListCell_;
  ASSERT(cellType_ == 1, "only implemented for cellType_ == 1");
  ASSERT(cellAtomCut_,
    "cellAtomCut must be on to use builNeighListCellAtomCut");
  const int iCellOfiAtom = atom2cell_[ipart];

  // add molecules in neighboring cells of iMol
  for (unsigned int i = 0; i < neighCell_[iCellOfiAtom].size(); ++i) {
    const int cell = neighCell_[iCellOfiAtom][i];
    for (unsigned int j = 0; j < cellList_[cell].size(); ++j) {
      neighListCell_.push_back(cellList_[cell][j]);
    }
  }
}

void Space::buildNeighListCell(const int iMol) {
  neighListCell_.clear();
  neighListChosen_ = &neighListCell_;
  ASSERT(cellType_ == 1, "only implemented for cellType_ == 1");
  ASSERT(!cellAtomCut_, "cellAtomCut must be off to use builNeighListCell");
  const int iCellOfiMol = imol2m(iMol);

  // add molecules in neighboring cells of iMol
  for (unsigned int i = 0; i < neighCell_[iCellOfiMol].size(); ++i) {
    const int cell = neighCell_[iCellOfiMol][i];
    for (unsigned int j = 0; j < cellList_[cell].size(); ++j) {
      neighListCell_.push_back(cellList_[cell][j]);
    }
  }
}

/**
 * turn off cell list
 */
void Space::cellOff() {
  cellType_ = 0;
}

void Space::eraseAtomFromCell_(const int ipart) {
  const int iCell = atom2cell_[ipart];
  int index = -1;
  vector<int> &icl = cellList_[iCell];
  for (unsigned int i = 0; i < icl.size(); ++i) {
    if (icl[i] == ipart) index = i;
  }
  icl.erase(icl.begin() + index);
  if (static_cast<int>(atom2cell_.size()) <= ipart) atom2cell_.resize(natom());
  atom2cell_[ipart] = -1;
}

void Space::eraseMolFromCell_(const int iMol) {
  const int iCell = mol2cell_[iMol];
  int index = -1;
  vector<int> &icl = cellList_[iCell];
  for (unsigned int i = 0; i < icl.size(); ++i) {
    if (icl[i] == iMol) index = i;
  }
  icl.erase(icl.begin() + index);
  if (static_cast<int>(mol2cell_.size()) <= iMol) mol2cell_.resize(nMol());
  mol2cell_[iMol] = -1;
}

void Space::addAtomtoCell_(const int ipart) {
  const int iCell = iatom2m(ipart);
  cellList_[iCell].push_back(ipart);
  if (static_cast<int>(atom2cell_.size()) <= ipart) {
    atom2cell_.resize(natom());
  }
  atom2cell_[ipart] = iCell;
}

void Space::addMoltoCell_(const int iMol) {
  const int iCell = imol2m(iMol);
  cellList_[iCell].push_back(iMol);
  if (static_cast<int>(mol2cell_.size()) <= iMol) mol2cell_.resize(nMol());
  // cout << "setting mol2cell_[" << iMol << "] = " << iCell << endl;
  mol2cell_[iMol] = iCell;
}

/**
 * updates cell for iMol, if, for example, it moves
 */
void Space::updateCellofiMol(const int iMol  //!< molecule to update
  ) {
  if (cellAtomCut_) {
    for (int ipart = mol2part_[iMol]; ipart < mol2part_[iMol+1]; ++ipart) {
      const int iCelln = iatom2m(ipart);
      const int iCello = atom2cell_[ipart];
      if (iCello != iCelln) {
        // cout << "moving from cells " << iCello << " -> " << iCelln << endl;
        eraseAtomFromCell_(ipart);
        addAtomtoCell_(ipart);
      }
    }
  } else {
    const int iCelln = imol2m(iMol);
    const int iCello = mol2cell_[iMol];
    if (iCello != iCelln) {
      // cout << "moving from cells " << iCello << " -> " << iCelln << endl;
      eraseMolFromCell_(iMol);
      addMoltoCell_(iMol);
    }
  }
}

/**
 * updates cell for iMol, if, for example, it moves
 */
void Space::updateCellofallMol() {
  if (cellAtomCut_) {
    for (int ipart = 0; ipart < natom(); ++ipart) {
      const int iCelln = iatom2m(ipart);
      const int iCello = atom2cell_[ipart];
      if (iCello != iCelln) {
        // cout << "moving from cells " << iCello << " -> " << iCelln << endl;
        eraseAtomFromCell_(ipart);
        addAtomtoCell_(ipart);
      }
    }
  } else {
    for (int iMol = 0; iMol < nMol(); ++iMol) {
      const int iCelln = imol2m(iMol);
      const int iCello = mol2cell_[iMol];
      if (iCello != iCelln) {
        // cout << "moving from cells " << iCello << " -> " << iCelln << endl;
        eraseMolFromCell_(iMol);
        addMoltoCell_(iMol);
      }
    }
  }
}

// stores current cell list, rebuilds, and compares
int Space::checkCellList() {
  int cellMatch = 1;

  // check mol2cell_
  if (cellAtomCut_) {
    ASSERT(static_cast<int>(atom2cell_.size()) == natom(),
      "atom2cell size(" << atom2cell_.size() << ") doesn't match number"
      << "of atoms(" << natom() << ")");
  } else {
    ASSERT(static_cast<int>(mol2cell_.size()) == nMol(),
    "mol2cell size(" << mol2cell_.size() << ") doesn't match number of"
    << "molecules(" << nMol() << ")");
    for (unsigned int i = 0; i < mol2cell_.size(); ++i) {
      ASSERT(mol2cell_[i] < nCell_, "mol2cell_[" << i << " has cell("
        << mol2cell_[i] << ") that doesn't exist");
      ASSERT(mol2cell_[i] >= 0, "mol2cell_[" << i << "] has -1 val");
    }
  }
  vector<int> mol2cell = mol2cell_;
  vector<int> atom2cell = atom2cell_;
  vector<vector<int> > cellList = cellList_;
  buildCellList();

  if (cellAtomCut_) {
    for (unsigned int i = 0; i < atom2cell_.size(); ++i) {
      ASSERT(atom2cell_[i] == atom2cell[i], "old atom2cell[" << i << "]="
        << atom2cell[i] << " doesn't match new atom2cell[" << i << "]="
        << atom2cell_[i] << " for particle with position " << x(i, 0)
        << " " << x(i, 1) << " " << x(i, 2));
    }
  } else {
    for (unsigned int i = 0; i < mol2cell_.size(); ++i) {
      ASSERT(mol2cell_[i] == mol2cell[i], "old mol2cell[" << i << "]="
        << mol2cell[i] << " doesn't match new mol2cell[" << i << "]="
        << mol2cell_[i] << "].");
    }
  }
  for (unsigned int i = 0; i < cellList_.size(); ++i) {
    // sort on-fly cell list for comparison
    std::sort(cellList[i].begin(), cellList[i].end());

    for (unsigned int j = 0; j < cellList_[i].size(); ++j) {
      ASSERT(cellList_[i][j] == cellList[i][j], "old cellList[" << i
        << "][" << j << "]=" << cellList[i][j]
        << " doesn't match new cellList[" << i << "][" << j << "]="
        << cellList_[i][j] << "].");
    }
  }
  return cellMatch;
}

void Space::initLMPData(const std::string fileName, const int nTypesExist) {
  // open LAMMPS data file
  std::ifstream file(fileName.c_str());
  ASSERT(file.good(), "cannot find lammps DATA file " << fileName);

  // skip all lines beginning with the character "#"
  skipCharsInFile('#', file);

  // read number of atoms
  string line, descript, descript2;
  int natoms = -1;
  file >> natoms >> descript;
  while (descript.compare("atoms") != 0) {
    file >> natoms >> descript;
  }

  // read next line, if it is number of bonds, then record.
  // Else, it is number of atom types
  int nBonds, nAngles, natype = 0;
  file >> nBonds >> descript;
  if (descript.compare("bonds") == 0) {
    bondList_.resize(nBonds);

    // read next line, if it is number of angles, then record.
    // Else, it is number of atom types
    file >> nAngles >> descript;
    if (descript.compare("angles") == 0) {
      angleList_.resize(nAngles);
    } else if (descript.compare("atom") == 0) {
      natype = nAngles;
    } else {
      ASSERT(0, "unrecognized lammps DATA format for file " << fileName);
    }
  } else if (descript.compare("atom") == 0) {
    natype = nBonds;
  } else {
    ASSERT(0, "unrecognized lammps DATA format for file " << fileName);
  }

  // read number of atom types, if not already done so
  if (natype == 0) {
    file >> natype >> descript >> descript2;
    while (descript.compare("atom") != 0) {
      file >> natype >> descript >> descript2;
    }
  } else {
    file >> descript2;
  }
  // cout << "extra " << descript << endl;
  nType_.resize(natype, 0);

  // cout << "nBonds " << bondList_.size() << " nA " << angleList_.size()
  //      << " nat " << natype << " nbt " << bondParam_.size() << " nat "
  //      << bondParam_.size() << endl;

  // read number of bond and angle types
  int nbondtypes, nangletypes;
  if (static_cast<int>(bondList_.size()) != 0) {
    file >> nbondtypes >> descript >> descript;
    bondParam_.resize(nbondtypes);
    // cout << "a " << nbondtypes << " " << descript  << endl;
  }
  if (static_cast<int>(angleList_.size()) != 0) {
    file >> nangletypes >> descript >> descript;
    angleParam_.resize(nangletypes);
  }

  // read bond parameters
  if (static_cast<int>(bondList_.size()) != 0) {
    readUntil("Bond Coeffs", file);
    double k, l0;
    int i;
    for (unsigned int j = 0; j < bondParam_.size(); ++j) {
      file >> i >> k >> l0;
      vector<double> b;
      b.push_back(k);
      b.push_back(l0);
      bondParam_[j] = b;
      getline(file, line);
    }
  }

  // read angle parameters
  if (static_cast<int>(angleList_.size()) != 0) {
    readUntil("Angle Coeffs", file);
    double k, t0;
    int i;
    for (unsigned int j = 0; j < angleParam_.size(); ++j) {
      file >> i >> k >> t0;
      vector<double> a;
      a.push_back(k);
      a.push_back(t0/180*PI);
      angleParam_[j] = a;
      getline(file, line);
    }
  }
  // cout << nBonds << " nA " << nAngles << " nat " << natype << " nbt "
  //      << bondParam_.size() << " nat " << bondParam_.size() << endl;

  // read until Atoms section
  readUntil("Atoms", file);

  // read Atoms section
  mol2part_.clear();
  listAtoms_.clear();
  int iatom, imol, itype, imolprev = 0;
  double q;
  vector<double> xtmp(2*dimen_);
  std::string cm, typetmp;
  for (int i = 0; i < natoms; ++i) {
    file >> iatom >> imol >> itype >> q;
    for (int dim = 0; dim < 2*dimen_; ++dim) {
      file >> xtmp[dim];
    }

    // add new atom
    for (int dim = 0; dim < dimen_; ++dim) {
      x_.push_back(xtmp[dim]);
    }
    mol_.push_back(imol - 1);
    type_.push_back(itype - 1 + nTypesExist);
    ++nType_[itype - 1];
    listAtoms_.push_back(i);

    // new molecule?
    if ( (i == 0) || (imol != imolprev) ) {
      imolprev = imol;

      // add new molecule
      mol2part_.push_back(iatom);
      moltype_.push_back(fileName);
      // moltype_.push_back(trim(".", typetmp.c_str()));
    }

    getline(file, line);
  }
  mol2part_.push_back(natom());
  xMolGen();
  if (natom() != 1) {
    qMolInit();
  }

  // read list of bonds
  if (static_cast<int>(bondList_.size()) != 0) {
    readUntil("Bonds", file);
    int iBond, iType, a1, a2;
    for (int i = 0; i < static_cast<int>(bondList_.size()); ++i) {
      file >> iBond >> iType >> a1 >> a2;
      // cout << "bond " << iType << " " << a1 << " " << a2 << endl;
      vector<int> b;
      b.push_back(iType-1);
      b.push_back(a1-1);
      b.push_back(a2-1);
      bondList_[i] = b;
    }
  }

  // read list of angles
  if (static_cast<int>(angleList_.size()) != 0) {
    readUntil("Angles", file);
    int iAngle, iType, a1, a2, a3;
    for (unsigned int i = 0; i < angleList_.size(); ++i) {
      file >> iAngle >> iType >> a1 >> a2 >> a3;
      vector<int> b;
      b.push_back(iType-1);
      b.push_back(a1-1);
      b.push_back(a2-1);
      b.push_back(a3-1);
      angleList_[i] = b;
    }
  }
}

int Space::checkSizes() {
  bool er = false;
  std::ostringstream ermesg;

  // per atom variables
  if (static_cast<int>(type_.size()) != natom()) {
    er = true;
    ermesg << "type(" << type_.size() << ") doesn't match natom("
           << natom() << ")." << endl;
  }
  if (static_cast<int>(mol_.size()) != natom()) {
    er = true;
    ermesg << "mol(" << mol_.size() << ") doesn't match natom(" << natom()
           << ")." << endl;
  }
  if (static_cast<int>(listAtoms_.size()) != natom()) {
    er = true;
    ermesg << "listAtoms(" << listAtoms_.size() << ") doesn't match natom("
           << natom() << ")." << endl;
  }

  // per molecule variables
  if (static_cast<int>(mol2part_.size()) - 1 != nMol()) {
    er = true;
    ermesg << "mol2part(" << mol2part_.size() - 1
           << ") doesn't match nMol(" << nMol() << ")." << endl;
  }
  if (static_cast<int>(moltype_.size()) != nMol()) {
    er = true;
    ermesg << "moltype(" << moltype_.size() << ") doesn't match nMol("
           << nMol() << ")." << endl;
  }
  if (static_cast<int>(listMols_.size()) != nMol()) {
    er = true;
    ermesg << "listMols(" << listMols_.size() << ") doesn't match nMol("
           << nMol() << ")." << endl;
  }
  if (sphereSymMol_ == false) {
    if (static_cast<int>(qMol_.size())/qdim_ != nMol()) {
      er = true;
      ermesg << "qMol(" << qMol_.size()/qdim_ << ") doesn't match nMol("
             << nMol() << ")." << endl;
    }
    if (static_cast<int>(xMolRef_.size()) != nMol()) {
      er = true;
      ermesg << "xMolRef(" << xMolRef_.size() << ") doesn't match nMol("
             << nMol() << ")." << endl;
    }
  } else {
    if (static_cast<int>(qMol_.size()) != 0) {
      er = true;
      ermesg << "qMol(" << qMol_.size() << ") should be zero for"
             << "spherically symmetric particle" << endl;
    }
    if (static_cast<int>(xMolRef_.size()) != 0) {
      er = true;
      ermesg << "xMolRef(" << xMolRef_.size() << ") should be zero"
             << "for spherically symmetric particle" << endl;
    }
  }
  if (cellType_ > 0) {
    if (cellAtomCut_) {
      if (static_cast<int>(atom2cell_.size()) != natom()) {
        er = true;
        ermesg << "atom2cell(" << atom2cell_.size() << ") doesn't match"
               << "natom(" << natom() << ")." << endl;
      }
      int natoms = 0;
      for (unsigned int i = 0; i < cellList_.size(); ++i) {
        for (unsigned int j = 0; j < cellList_[i].size(); ++j) {
          ++natoms;
        }
      }
      if (natoms != natom()) {
        er = true;
        ermesg << "cellList(" << natoms << ") doesn't match natom("
               << natom() << ")." << endl;
      }
    } else {
      if (static_cast<int>(mol2cell_.size()) != nMol()) {
        er = true;
        ermesg << "mol2cell(" << mol2cell_.size() << ") doesn't match nMol("
               << nMol() << ")." << endl;
      }
      int natoms = 0;
      for (unsigned int i = 0; i < cellList_.size(); ++i) {
        for (unsigned int j = 0; j < cellList_[i].size(); ++j) {
          ++natoms;
        }
      }
      if (natoms != nMol()) {
        er = true;
        ermesg << "cellList(" << natoms << ") doesn't match nMol("
               << nMol() << ")." << endl;
      }
    }
    if (cellType_ == 1) {
      for (unsigned int i = 0; i < neighCell_.size(); ++i) {
        if (neighCell_[i].size() != pow(3, dimen_)) {
          er = true;
          ermesg << "For cellType(" << cellType_ << ") there should be"
            << "3^dim-1(" << pow(3, dimen_) << ") neighbors, however,"
            << "neighCell_[" << i << "].size = " << neighCell_[i].size()
            << endl;
        }
      }
    }
  }

  if (er) {
    ASSERT(0, "size check failure" << endl << ermesg.str());
    return 0;
  }
  return 1;
}

int Space::mvec2m3d_(const int &i, const int &j, const int &k) const {
  const int mx = nCellVec_[0], my = nCellVec_[1], mz = nCellVec_[2];
  return (i+mx)%mx + mx*((j+my)%my + my*((k+mz)%mz));
}

int Space::mvec2m2d_(const int &i, const int &j) const {
  const int mx = nCellVec_[0], my = nCellVec_[1];
  return (i+mx)%mx + mx*((j+my)%my);
}

shared_ptr<Space> Space::findAddMolInList(const string typeStr) {
  bool match = false;
  for (unsigned int iaml = 0; iaml < addMolList_.size(); ++iaml) {
    if (addMolListType_[iaml].compare(typeStr) == 0) {
      if (match) {
        ASSERT(0, "double match in addMol");
      } else {
        match = true;
      }
      return addMolList_[iaml];
    }
  }
  ASSERT(match, "unrecognized molecule type(" << typeStr
         << ") for addMol function of space.cc");
  return shared_ptr<Space>();
}

int Space::findAddMolListIndex(const string typeStr) {
  bool match = false;
  for (unsigned int iaml = 0; iaml < addMolList_.size(); ++iaml) {
    if (addMolListType_[iaml].compare(typeStr) == 0) {
      return iaml;
    }
  }
  ASSERT(match, "unrecognized molecule type(" << typeStr
         << ") for addMol function of space.cc");
  return -1;
}

void Space::scaleMol(const int iMol, const vector<double> bondLengths) {
  vector<vector<double> > xMolRef = xMolRef_[iMol];

  // for each atom in molecule
  for (unsigned int iAtom = 0; iAtom < xMolRef.size(); ++iAtom) {
    double rinv = 1., r2 = 0.;
    for (int dim = 0; dim < dimen_; ++dim) {
      r2 += pow(xMolRef[iAtom][dim] - xMolRef[0][dim], 2);
    }
    if (r2 != 0.) rinv = 1./sqrt(r2);
    for (int dim = 0; dim < dimen_; ++dim) {
      xMolRef_[iMol][iAtom][dim] *= bondLengths[iAtom] * rinv;
    }
  }
  quat2pos(iMol);
}

vector<int> Space::tag2mpart() {
  vector<int> mpart;
  const int ipart = tag_[0];
  const int iMol = mol_[ipart];
  ASSERT(ipart == mol2part_[iMol], "ipart(" << ipart
         << ") != mol2part_[iMol=" << iMol << "]=" << mol2part_[iMol]);
  ASSERT(tag_.size() == 1, "tag size " << tag_.size());
  for (int i = ipart; i < mol2part_[iMol+1]; ++i) {
    mpart.push_back(i);
  }
  return mpart;
}

void Space::writeRestart(const char* fileName) {
  fileBackUp(fileName);
  std::ofstream file(fileName);

  // print spatial parameters
  file << "# id " << id_ << endl;
  file << "# dimen " << dimen_ << endl;
  for (int dim = 0; dim < dimen_; ++dim) {
    file << std::setprecision(std::numeric_limits<double>::digits10+2)
         << "# l" << dim << " " << l_[dim] << endl;
  }
  if (fabs(xyTilt_) > doubleTolerance) file << "# xyTilt " << xyTilt_ << endl;
  if (fabs(xzTilt_) > doubleTolerance) file << "# xzTilt " << xzTilt_ << endl;
  if (fabs(yzTilt_) > doubleTolerance) file << "# yzTilt " << yzTilt_ << endl;
  if (floppyBox_ != 0) file << "# floppyBoxFlag " << floppyBox_;
  if (maxlFlag_ != 0) {
    file << "# maxlFlag " << maxlFlag_;
    for (int dim = 0; dim < dimen_; ++dim) {
      file << "# maxboxl" << dim << " " << maxl_[dim] << endl;
    }
  }

  file << "# cellAtomCut " << cellAtomCut_ << endl;
  if (cellType_ > 0) {
    file << "# cellType " << cellType_ << endl;
    file << std::setprecision(std::numeric_limits<double>::digits10+2)
         << "# dCellMin " << dCellMin_ << endl;
  }

  // print addmolinits
  file << "# naddmolinits " << addMolListType_.size() << endl;
  for (unsigned int i = 0; i < addMolListType_.size(); ++i) {
    file << "# addmolinittype" << i << " " << addMolListType_[i] << endl;
  }

  // print type and number of each type of molecule
  int iMol = 0, iMolPrev = 0, nMolTypes = 0;
  while (iMol < nMol()) {
    string mtype = moltype_[iMol], mtypenext = mtype;
    // count number of consecutive molecules of each type
    while ( (iMol < nMol()) && (mtype.compare(mtypenext) == 0) ) {
      ++iMol;
      if (iMol < nMol()) mtypenext = moltype_[iMol];
    }
    file << "# moltype" << nMolTypes << " " << mtype << endl;
    file << "# nmoloft" << nMolTypes << " " << iMol-iMolPrev << endl;
    iMolPrev = iMol;
    ++nMolTypes;
  }
  file << "# nMolTypes " << nMolTypes << endl;

  // print something about tagged molecule?
  if (static_cast<int>(tag_.size()) > 0) {
    file << "# ntags " << tag_.size() << endl;
    for (unsigned int i = 0; i < tag_.size(); ++i) {
      file << "# tag" << i << " " << tag_[i] << endl;
    }
  }

  // print something about cluster analysis
  file << "# nClusterTypes " << clusterType_.size() << endl;
  for (unsigned int i = 0; i < clusterType_.size(); ++i) {
    file << "# clusterType" << i << " " << clusterType_[i] << endl;
  }

  if (equiMolar_ != 0) file << "# equiMolarFlag " << equiMolar_ << endl;

  // write random number generator state
  writeRngRestart(fileName);

  if (eulerFlag_ != 0) file << "# eulerFlag " << eulerFlag_ << endl;

  // print configuration
  if (sphereSymMol_) {
    // for each atom, print coordinates
    for (int iAtom = 0; iAtom < natom(); ++iAtom) {
      for (int dim = 0; dim < dimen_; ++dim) {
        file << std::setprecision(std::numeric_limits<double>::digits10+2)
             << x(iAtom, dim) << " ";
      }
      file << endl;
    }
  } else {
    // for each molecule, print position of first(pivot) atom, xMolRef, and qMol
    for (int iMol = 0; iMol < nMol(); ++iMol) {
      const int iAtom = mol2part_[iMol];
      for (int dim = 0; dim < dimen_; ++dim) {
        file << std::setprecision(std::numeric_limits<double>::digits10+2)
             << x(iAtom, dim) << " ";
      }
      file << endl;
      for (unsigned int ipart = 1; ipart < xMolRef_[iMol].size(); ++ipart) {
        for (int dim = 0; dim < dimen_; ++dim) {
          file << std::setprecision(std::numeric_limits<double>::digits10+2)
               << xMolRef_[iMol][ipart][dim] << " ";
        }
        file << endl;
      }
      for (int qdim = 0; qdim < qdim_; ++qdim) {
        file << std::setprecision(std::numeric_limits<double>::digits10+2)
             << qMol(iMol, qdim) << " ";
      }
      file << endl;
    }
  }
}

/**
 * flood fill algorithm to identify clusters based on atomic distance cutoff
 */
void Space::floodFill3d_(
  const int clusterNode,  //!< atom on cluster edge to grow
  const int clusterID,    //!< id of current cluster
  const double rCut       //!< distance based cutoff to define cluster
  ) {
  // declare variables for optimization
  const double xi = xcluster_[dimen_*clusterNode];
  const double yi = xcluster_[dimen_*clusterNode+1];
  const double zi = xcluster_[dimen_*clusterNode+2];
  const double rCut2 = rCut * rCut;
  double r2, dx, dy, dz, dx2, dy2, dz2;
  const double lx = l_[0], ly = l_[1], lz = l_[2], halflx = lx/2.,
    halfly = ly/2., halflz = lz/2.;

  // find all nearest neighbors not already assigned a cluster type
  // within rCut of clusterNode.
  // Ror each neighbor, recursively active another floodFill
  for (int i = 0; i < natom(); ++i) {
    if (cluster_[i] == -natom()) {
      // separation distance with periodic boundary conditions
      dx = xi - xcluster_[dimen_*i];
      dy = yi - xcluster_[dimen_*i+1];
      dz = zi - xcluster_[dimen_*i+2];
      dx2 = dy2 = dz2 = 0.;
      if (dx >  halflx) {
        dx -= lx;
        dx2 = -lx;
      }
      if (dx < -halflx) {
        dx += lx;
        dx2 = lx;
      }
      if (dy >  halfly) {
        dy -= ly;
        dy2 = -ly;
      }
      if (dy < -halfly) {
        dy += ly;
        dy2 = ly;
      }
      if (dz >  halflz) {
        dz -= lz;
        dz2 = -lz;
      }
      if (dz < -halflz) {
        dz += lz;
        dz2 = lz;
      }
      r2 = dx*dx + dy*dy + dz*dz;

      if (r2 < rCut2) {
        cluster_[i] = clusterID;
        const int imol = mol_[i];
        clusterMol_[imol] = clusterID;

        // move molecule in xcluster_ if pbc was used
        if ((dx2 != 0) || (dy2 != 0) || (dz2 != 0)) {
          for (int ipart = mol2part_[imol]; ipart < mol2part_[imol+1];
               ++ipart) {
            xcluster_[dimen_*ipart    ] -= dx2;
            xcluster_[dimen_*ipart + 1] -= dy2;
            xcluster_[dimen_*ipart + 2] -= dz2;
          }
        }

        floodFill3d_(i, clusterID, rCut);
      }
    }
  }
}

/**
 * flood fill algorithm to identify clusters based on atomic distance cutoff
 */
void Space::floodFill2d_(
  const int clusterNode,  //!< atom on cluster edge to grow
  const int clusterID,    //!< id of current cluster
  const double rCut       //!< distance based cutoff to define cluster
  ) {
  // declare variables for optimization
  const double xi = xcluster_[dimen_*clusterNode];
  const double yi = xcluster_[dimen_*clusterNode+1];
  const double rCut2 = rCut * rCut;
  double r2, dx, dy, dx2, dy2;
  const double lx = l_[0], ly = l_[1], halflx = lx/2., halfly = ly/2.;

  // find all nearest neighbors not already assigned a cluster type
  // within rCut of clusterNode
  // For each neighbor, recursively active another floodFill
  for (int i = 0; i < natom(); ++i) {
    if (cluster_[i] == -natom()) {
      // separation distance with periodic boundary conditions
      dx = xi - xcluster_[dimen_*i];
      dy = yi - xcluster_[dimen_*i+1];
      dx2 = dy2 = 0.;
      if (dx >  halflx) {
        dx -= lx;
        dx2 = -lx;
      }
      if (dx < -halflx) {
        dx += lx;
        dx2 = lx;
      }
      if (dy >  halfly) {
        dy -= ly;
        dy2 = -ly;
      }
      if (dy < -halfly) {
        dy += ly;
        dy2 = ly;
      }
      r2 = dx*dx + dy*dy;

      if (r2 < rCut2) {
        cluster_[i] = clusterID;
        const int imol = mol_[i];
        clusterMol_[imol] = clusterID;

        // move molecule in xcluster_ if pbc was used
        if ((dx2 != 0) || (dy2 != 0)) {
          for (int ipart = mol2part_[imol]; ipart < mol2part_[imol+1];
               ++ipart) {
            xcluster_[dimen_*ipart    ] -= dx2;
            xcluster_[dimen_*ipart + 1] -= dy2;
          }
        }

        floodFill2d_(i, clusterID, rCut);
      }
    }
  }
}

/**
 * flood fill algorithm to identify clusters based on atomic distance cutoff
 *  uses cell list
 * HWH CLEANUP: Alt? too much copy and paste
 */
void Space::floodFillCell3d_(
  const int clusterNode,  //!< atom on cluster edge to grow
  const int clusterID,    //!< id of current cluster
  const double rCut       //!< distance based cutoff to define cluster
  ) {
  // declare variables for optimization
  const double xi = xcluster_[dimen_*clusterNode];
  const double yi = xcluster_[dimen_*clusterNode+1];
  const double zi = xcluster_[dimen_*clusterNode+2];
  const double rCut2 = rCut * rCut;
  double r2, dx, dy, dz, dx2, dy2, dz2;
  const double lx = l_[0], ly = l_[1], lz = l_[2],
    halflx = lx/2., halfly = ly/2., halflz = lz/2.;
  const int iCell = atom2cell_[clusterNode];

  // find all nearest neighbors not already assigned a cluster type
  // within rCut of clusterNode
  // For each neighbor, recursively active another floodFillCell
  for (unsigned int nc = 0; nc < neighCell_[iCell].size(); ++nc) {
    const int jCell = neighCell_[iCell][nc];
    for (unsigned int j = 0; j < cellList_[jCell].size(); ++j) {
      const int i = cellList_[jCell][j];
      if (cluster_[i] == -natom()) {
        // separation distance with periodic boundary conditions
        dx = xi - xcluster_[dimen_*i];
        dy = yi - xcluster_[dimen_*i+1];
        dz = zi - xcluster_[dimen_*i+2];
        dx2 = dy2 = dz2 = 0.;
        if (dx >  halflx) {
          dx -= lx;
          dx2 = -lx;
        }
        if (dx < -halflx) {
          dx += lx;
          dx2 = lx;
        }
        if (dy >  halfly) {
          dy -= ly;
          dy2 = -ly;
        }
        if (dy < -halfly) {
          dy += ly;
          dy2 = ly;
        }
        if (dz >  halflz) {
          dz -= lz;
          dz2 = -lz;
        }
        if (dz < -halflz) {
          dz += lz;
          dz2 = lz;
        }
        r2 = dx*dx + dy*dy + dz*dz;

        if (r2 < rCut2) {
          cluster_[i] = clusterID;
          const int imol = mol_[i];
          clusterMol_[imol] = clusterID;

          // move molecule in xcluster_ if pbc was used
          if ((dx2 != 0) || (dy2 != 0) || (dz2 != 0)) {
            for (int ipart = mol2part_[imol]; ipart < mol2part_[imol+1];
                 ++ipart) {
              xcluster_[dimen_*ipart    ] -= dx2;
              xcluster_[dimen_*ipart + 1] -= dy2;
              xcluster_[dimen_*ipart + 2] -= dz2;
            }
          }

          floodFillCell3d_(i, clusterID, rCut);
        }
      }
    }
  }
}

/**
 * prefil cluster vars
 */
void Space::prefilClusterVars_() {
  // prefill cluster vector with -natom() if included and -natom-1 if excluded
  cluster_.resize(natom());
  for (int i = 0; i < natom(); ++i) {
    if (findInList(type_[i], clusterType_)) {
      cluster_[i] = -natom();
    } else {
      cluster_[i] = -natom()-1;
    }
  }
  clusterMol_.resize(nMol());
  std::fill(clusterMol_.begin(), clusterMol_.end(), -1);
  xcluster_ = x_;
}

void Space::updateClusters(const double rCut) {
  prefilClusterVars_();
  int nClusters = 0;
  for (int i = 0; i < natom(); ++i) {
    if (cluster_[i] == -natom()) {
      if ( (cellType_ == 1) && (rCut <= dCellMin()) && (dimen_ == 3) ) {
        floodFillCell3d_(i, nClusters, rCut);
      } else {
        if (dimen_ == 3) {
          floodFill3d_(i, nClusters, rCut);
        } else if (dimen_ == 2) {
          floodFill2d_(i, nClusters, rCut);
        } else {
          ASSERT(0, "floodFill algorithm not implemented for dimen:" << dimen_);
        }
      }
      ++nClusters;
    }
  }

  // generate cluster atom list and cluster sizes
  updateClusterVars(nClusters);
}

void Space::updateClusterVars(const int nClusters) {
  ASSERT(nClusters != 0, "no clusters found. Did you use addTypeForCluster"
    << "function to define clusters? Or is natom(" << natom() << ") == 0?");
  clusterSizes_.resize(nClusters);
  std::fill(clusterSizes_.begin(), clusterSizes_.end(), 0);
  clusterList_.clear();
  clusterList_.resize(nClusters);
  for (int i = 0; i < natom(); ++i) {
    ASSERT(cluster_[i] != -natom(), "cluster is -natom");
    if (cluster_[i] != -natom()-1) {
      ++clusterSizes_[cluster_[i]];
    }
    ASSERT(clusterMol_[mol_[i]] != -1,
      "clusterMol[mol[" << i << "]=" << mol_[i] << "]=-1");
    clusterList_[clusterMol_[mol_[i]]].push_back(i);
  }

  // accumulate cluster size and number statistics
  clusterSizeAccVec_.accumulate(nMol(), clusterAvSize());
  clusterNumAccVec_.accumulate(nMol(), nClusters);
  AccumulatorVec csdist;
  int nfree = 0;
  for (int i = 0; i < nClusters; ++i) {
    csdist.accumulate(clusterSizes_[i], 1);
    if (clusterSizes_[i] <= preMicellarAgg_) ++nfree;
  }
  for (int i = 0; i < csdist.size(); ++i) {
    clusterSizeDistribution_.accumulate(i, csdist.vec(i).sum());
    if (peStore_ != -1) {
      const double csdvs = csdist.vec(i).sum();
      clusterSizeDistributionU_.accumulate(i, csdvs*peStore_);
      clusterSizeDistributionU2_.accumulate(i, csdvs*peStore_*peStore_);
    }
  }
  freeMon_.accumulate(nfree/vol());
}

void Space::contact2cluster(
  vector<vector<int> > contact,
  vector<vector<vector<double> > > contactpbc
  ) {
  ASSERT(0, "xcluster issue. use contact2clusterAlt");
  prefilClusterVars_();
  int nClusters = 0;
  for (int i = 0; i < natom(); ++i) {
    if (cluster_[i] == -natom()) {
      cluster_[i] = nClusters;
      clusterMol_[mol_[i]] = nClusters;
      floodFillContact_(i, nClusters, &contact, &contactpbc);
      ++nClusters;
    }
  }

  updateClusterVars(nClusters);
}

void Space::floodFillContact_(const int clusterNode,
  const int clusterID,
  vector<vector<int> > *contactPtr,
  vector<vector<vector<double> > > *contactpbcPtr) {
  vector<vector<int> >& contact = *contactPtr;
  vector<vector<vector<double> > >& contactpbc = *contactpbcPtr;
  const int iMol = mol_[clusterNode];

  for (int i = 0; i < natom(); ++i) {
    if (cluster_[i] == -natom()) {
      const int jMol = mol_[i];
      if (contact[iMol][jMol] == 1) {
        cluster_[i] = clusterID;
        clusterMol_[jMol] = clusterID;
        for (int ipart = mol2part_[jMol]; ipart < mol2part_[jMol+1]; ++ipart) {
          for (int dim = 0; dim < dimen_; ++dim) {
            xcluster_[dimen_*ipart+dim] -= contactpbc[iMol][jMol][dim];
          }
        }
        floodFillContact_(i, clusterID, contactPtr, contactpbcPtr);
      }
    }
  }
}

void Space::contact2clusterAlt(
  vector<vector<int> > contact,
  vector<vector<vector<double> > > contactpbc
  ) {
  prefilClusterVars_();
  int nClusters = 0;
  percolation_ = 0;
  for (int i = 0; i < natom(); ++i) {
    if (cluster_[i] == -natom()) {
      cluster_[i] = nClusters;
      clusterMol_[mol_[i]] = nClusters;
      vector<vector<int> > image(natom(), vector<int>(dimen_, 0));
      floodFillContactAlt_(i, nClusters, &contact, &contactpbc, &image);
      ++nClusters;
    }
  }

  updateClusterVars(nClusters);

//  // print xcluster
//  cout << "xclus" << endl;
//  for (int iatom = 0; iatom < natom(); ++iatom) {
//    for (int dim = 0; dim < dimen_; ++dim) {
//      cout << xcluster_[dimen_*iatom+dim] << " ";
//    }
//    cout << endl;
//  }
}

/**
 * flood fill algorithm with contact map
 * HWH CLEANUP: Alt? too much copy and paste
 */
void Space::floodFillContactAlt_(const int clusterNode,
  const int clusterID,
  vector<vector<int> > *contactPtr,
  vector<vector<vector<double> > > *contactpbcPtr,
  vector<vector<int> > *image) {
  vector<vector<int> >& contact = *contactPtr;
  vector<vector<vector<double> > >& contactpbc = *contactpbcPtr;
  const int iMol = mol_[clusterNode];

  for (int i = 0; i < natom(); ++i) {
    // if cluster==-natom, no cluster found yet, so check for contact
    // or, if == clusterID (already found), check image for percolation
    if ( (cluster_[i] == clusterID) || (cluster_[i] == -natom()) ) {
      const int jMol = mol_[i];
      int index = -1;
      if (findInList(jMol, contact[iMol], index)) {
        // compute the current image
        vector<int> currentImage(dimen_);
//        cout << "dpbc ";
//        int jindex = -1;
//        ASSERT(findInList(iMol, contact[jMol], jindex), "iMol/jMol reciprocity");
        for (int dim = 0; dim < dimen_; ++dim) {
          const double dpbc = contactpbc[iMol][index][dim];
          if (fabs(dpbc) < DTOL) {
            currentImage[dim] = (*image)[clusterNode][dim];
          } else if (dpbc > 0) {
            currentImage[dim] = (*image)[clusterNode][dim] - 1;
          } else if (dpbc < 0) {
            currentImage[dim] = (*image)[clusterNode][dim] + 1;
          }
//          cout << dpbc << " ";
        }
//        cout << endl;
//        cout << "Current image: " << vec2str(currentImage) << endl;

        // if in contact but not listed, add (i,jMol) to cluster
        if (cluster_[i] == -natom()) {

          // cout << "pbc " << vec2str(contactpbc[iMol][index]) << endl;
          // cout << "Current part: i " << i << " mol " << jMol << " node " << clusterNode << " mol " << iMol << " x " << x(i, 2) << endl;
          cluster_[i] = clusterID;
          clusterMol_[jMol] = clusterID;
          for (int ipart = mol2part_[jMol]; ipart < mol2part_[jMol+1]; ++ipart) {
            for (int dim = 0; dim < dimen_; ++dim) {
              xcluster_[dimen_*ipart+dim] += currentImage[dim]*l_[dim];
              //xcluster_[dimen_*ipart+dim] -= contactpbc[iMol][index][dim];
            }
          }

          // store the current image
          for (int dim = 0; dim < dimen_; ++dim) {
            (*image)[i][dim] = currentImage[dim];
          }

          floodFillContactAlt_(i, clusterID, contactPtr, contactpbcPtr, image);

        // if contact already found previously, check image for percolation
        } else {
//          cout << "already found" << endl;
//          cout << "Current part: i " << i << " mol " << jMol << " node " << clusterNode << " mol " << iMol << " x " << x(i, 2) << endl;
//          cout << "Previous image: " << vec2str((*image)[i]) << endl;
          for (int dim = 0; dim < dimen_; ++dim) {
            if ((*image)[i][dim] != currentImage[dim]) {
//              cout << "**" << endl << "PERCOLATION found!" << endl;
              //cout << "node " << clusterNode << endl;
              //cout << "i " << i << endl;
              percolation_ = 1;
            }
          }
        }
      }
    }
  }
}

void Space::delTypePart(const int type) {
  int i = 0;
  while (i != natom()) {
    if (type_[i] == type) {
      delPart(i);
    } else {
      ++i;
    }
  }
}

void Space::swapPositions(Space *space) {
  ASSERT(natom() == space->natom(), "natom(" << natom() << ") of space id "
    << id_ << " doesn't match natom(" << space->natom() << ") of space id "
    << space->id() << ")");
  ASSERT(dimen_ == space->dimen(), "dimen(" << dimen() << ") of space id "
    << id_ << " doesn't match dimen(" << space->dimen() << ") of space id "
    << space->id() << ")");
  vector<double> x = x_;
  for (unsigned int i = 0; i < x.size(); ++i) {
    double xtmp = x_[i];
    x_[i] = space->x_[i];
    space->x_[i] = xtmp;
  }
  vector<double> qMol = qMol_;
  for (unsigned int i = 0; i < qMol.size(); ++i) {
    double qMoltmp = qMol_[i];
    qMol_[i] = space->qMol_[i];
    space->qMol_[i] = qMoltmp;
  }
  vector<vector<vector<double> > > xMolRef = xMolRef_;
  for (unsigned int i = 0; i < xMolRef.size(); ++i) {
    for (unsigned int j = 0; j < xMolRef[i].size(); ++j) {
      for (unsigned int k = 0; k < xMolRef[i][j].size(); ++k) {
        double xMolReftmp = xMolRef_[i][j][k];
        xMolRef_[i][j][k] = space->xMolRef_[i][j][k];
        space->xMolRef_[i][j][k] = xMolReftmp;
      }
    }
  }
}

double Space::maxMolDist() {
  double max = 0.;
  xMolGen();

  // loop through all existing molecules
  for (unsigned int i = 0; i < xMol_.size(); ++i) {
    for (unsigned int j = 1; j < xMol_[i].size(); ++j) {
      const double r2 = rsq(xMol_[i][0], xMol_[i][j]);
      if (r2 > max) max = r2;
    }
  }

  // loop through potential molecules to add
  for (unsigned int i = 0; i < addMolList_.size(); ++i) {
    const double r2 = pow(addMolList_[i]->maxMolDist(), 2.);
    if (r2 > max) max = r2;
  }
  return sqrt(max);
}

/**
 * initialize cell list for atom cutoff
 */
void Space::initCellAtomCut(const int flag) {
  if (flag == 0) {
    cellAtomCut_ = false;
    neighListChosen_ = &listMols_;
  } else {
    cellAtomCut_ = true;
    neighListChosen_ = &listAtoms_;
  }
}

vector<int> Space::imol2mpart(const int iMol) {
  vector<int> mpart;
  ASSERT(iMol < nMol(),
    "in imol2mpart, iMol(" << iMol << ") >= nMol(" << nMol());
  for (int ipart = mol2part_[iMol]; ipart < mol2part_[iMol+1]; ++ipart) {
    mpart.push_back(ipart);
  }
  return mpart;
}

/**
 * Attempt translation of molecule
 */
void Space::transMol(const int iMol,    //!< molecule to translate
  const vector<double> &r   //!< translation vector
  ) {
  const vector<int> mpart = imol2mpart(iMol);
  ASSERT(static_cast<int>(r.size()) == dimen_,
    "size of r(" << r.size() << ") != dimen(" << dimen_ << ")");

  for (unsigned int i = 0; i < mpart.size(); ++i) {
    const int ipart = mpart[i];
    for (int dim = 0; dim < dimen_; ++dim) {
      x_[dimen_*ipart+dim] += r[dim];
    }
  }

  // wrap inside box
  wrap(mpart);
}

vector<vector<double> > Space::inertialTensor(const vector<int> mpart) {
  ASSERT(dimen_ == 3,
    "dimen(" << dimen_ << ") must be 3 in inertialTensor computation");
  vector<vector<double> > tensor(dimen_, vector<double>(dimen_, 0.));
  const int natom = static_cast<int>(mpart.size());
  vector<double> rcm = rcom(mpart);
  for (int i = 0; i < natom; ++i) {
    const int ipart = mpart[i];
    const double dx = x_[dimen_*ipart  ] - rcm[0];
    const double dy = x_[dimen_*ipart+1] - rcm[1];
    const double dz = x_[dimen_*ipart+2] - rcm[2];
    tensor[0][0] += dy*dy + dz*dz;
    tensor[1][1] += dx*dx + dz*dz;
    tensor[2][2] += dx*dx + dy*dy;
    tensor[0][1] = tensor[1][0] += dx*dy;
    tensor[1][2] = tensor[2][1] += dy*dz;
    tensor[0][2] = tensor[2][0] += dx*dz;
  }
  return tensor;
//    // check that moment of inertia is diagonal
}

vector<double> Space::rcom(const vector<int> mpart) {
  vector<double> rcom(dimen_, 0.);
  const int natom = static_cast<int>(mpart.size());
  for (int i = 0; i < natom; ++i) {
    for (int dim = 0; dim < dimen_; ++dim) {
      rcom[dim] += x_[dimen_*mpart[i]+dim] / static_cast<double>(natom);
    }
  }
  return rcom;
}

vector<int> Space::mpart2mmol(const vector<int> mpart) {
  vector<int> molList;
  int iMolPrev = -1;
  for (unsigned int ipart = 0; ipart < mpart.size(); ++ipart) {
    const int iMol = mol_[mpart[ipart]];
    if (iMol != iMolPrev) {
      molList.push_back(iMol);
    }
    iMolPrev = iMol;
  }
  return molList;
}

void Space::printClusterStat(const char* fileName) {
  // print size distribution
  stringstream ss;
  ss << fileName << "sizedist";
  std::ofstream outFile2(ss.str().c_str());
  outFile2 << "# clusterSize probability" << endl;
  vector<double> p(clusterSizeDistribution_.size());
  const double sum = clusterSizeDistribution_.sum();
  for (unsigned int i = 0; i < p.size(); ++i) {
    outFile2 << i << " " << clusterSizeDistribution_.vec(i).sum() / sum << endl;
  }

  // find first minimum in size distribution to free monomer
  // and premicellar concentration (cmc)
  vector<int> min = findLocalMinima(p, 1);

  // print statistics as function of number of molecules
  std::ofstream outFile(fileName);
  outFile << "# natom avClusterNum avClusterSize ";
  double nfree = 0;
  if (min.size() > 1) {
    outFile << "freemon" << endl;
    for (int i = 0; i <= min[1]; ++i) {
      nfree += clusterSizeDistribution_.vec(i).average();
    }
  }
  for (int i = 0; i < clusterSizeAccVec_.size(); ++i) {
    if (clusterNumAccVec_.vec(i).nValues() > 0) {
      outFile << i << " " << clusterNumAccVec_.vec(i).average() << " "
              << clusterSizeAccVec_.vec(i).average() << " ";
      if (static_cast<int>(min.size()) != 0) {
        outFile << nfree/vol();
      }
      outFile << endl;
    }
  }
}

void Space::xClusterGen() {
  xcluster_ = x_;
  for (int ic = 0; ic < nClusters(); ++ic) {
    const int iMol = mol_[clusterList_[ic][0]];
    const int ipart = mol2part_[iMol];
    for (int jMol = 0; jMol < nMol(); ++jMol) {
      if ( (jMol != iMol) && (clusterMol_[jMol] == ic) ) {
        const int jpart = mol2part_[jMol];
        vector<double> rij(dimen_);
        for (int dim = 0; dim < dimen(); ++dim) {
          rij[dim] = x(ipart, dim) - x(jpart, dim);
        }
        const vector<double> p = pbc(rij);
        for (int j = jpart; j < mol2part_[jMol+1]; ++j) {
          for (int dim = 0; dim < dimen(); ++dim) {
            xcluster_[dimen_*j+dim] -= p[dim];
          }
        }
      }
    }
  }
}

void Space::xClusterShape() {
  clusterAsphericity_.clear();
  clusterAcylindricity_.clear();
  clusterRelShapeAniso_.clear();
  clusterRg_.clear();

  // loop through each cluster
  for (int ic = 0; ic < nClusters(); ++ic) {
    const int cSize = static_cast<int>(clusterList_[ic].size());

    // find the COM of cluster
    vector<double> xcCOM(dimen_, 0.);
    for (int j = 0; j < cSize; ++j) {
      const int ipart = clusterList_[ic][j];
      for (int dim = 0; dim < dimen_; ++dim) {
        xcCOM[dim] += xcluster_[dimen_*ipart + dim] /
          static_cast<double>(cSize);
      }
    }

    // compute the gyration tensor of cluster
    vector<vector<double> > xcGy(dimen_, vector<double>(dimen_, 0.));
    for (int j = 0; j < cSize; ++j) {
      const int ipart = clusterList_[ic][j];
      for (int idim = 0; idim < dimen_; ++idim) {
        const double xi = xcluster_[dimen_*ipart + idim] - xcCOM[idim];
        for (int jdim = 0; jdim < dimen_; ++jdim) {
          const double xj = xcluster_[dimen_*ipart + jdim] - xcCOM[jdim];
          xcGy[idim][jdim] += xi*xj / static_cast<double>(cSize);
        }
      }
    }


    // find eigenvalues of gyration tensor
    vector<double> evalues(dimen_);
    vector<vector<double> > evectors(dimen_, vector<double>(dimen_));
    jacobi(xcGy, evalues, evectors);

    // compute the shape parameters
    ASSERT(dimen_ == 3,
      "asphericity computation assumes dim(" << dimen_ << ") = 3");
    std::sort(evalues.begin(), evalues.end());
    const double lx2 = evalues[0];
    const double ly2 = evalues[1];
    const double lz2 = evalues[2];
    const double rg2 = lx2+ly2+lz2;
    const double asphr = lz2 - 0.5*(lx2 + ly2);
    const double acyln = ly2 - lx2;
    const double relsh = (asphr*asphr + (3./4.)*acyln*acyln)/(rg2*rg2);
    clusterAsphericity_.push_back(asphr);
    clusterAcylindricity_.push_back(acyln);
    clusterRelShapeAniso_.push_back(relsh);
    clusterRg_.push_back(pow(rg2, 0.5));
  }
}

#ifdef XDRFILE_H_
int Space::readXTC(const char* fileName,
  XDRFILE* trjFileXDR) {
  int endXTC = 0;
  string fileNameStr(fileName);

  // check atoms
  int natoms_xtc;
  char * fn_xtc = new char[fileNameStr.size() + 1];
  std::copy(fileNameStr.begin(), fileNameStr.end(), fn_xtc);
  fn_xtc[fileNameStr.size()] = '\0';
  int result_xtc = read_xtc_natoms(fn_xtc, &natoms_xtc);
  ASSERT(exdrOK == result_xtc,
    "cannot read natoms, " << natoms_xtc << " from xtc file: " << fileNameStr);
  ASSERT(natom() == natoms_xtc, "number of atoms(" << natom()
    << ") does not match number of atoms in xtc, " << natoms_xtc
    << " from xtc file: " << fileNameStr);

  // set coordinates and box length
  int step_xtc;
  float time_xtc;
  matrix box_xtc;
  rvec *x_xtc;
  x_xtc = reinterpret_cast<rvec *>(calloc(natoms_xtc, sizeof(x_xtc[0])));
  // x_xtc = (rvec *)calloc(natoms_xtc, sizeof(x_xtc[0]));
  float prec_xtc = 1000.0;
  result_xtc = read_xtc(trjFileXDR, natoms_xtc, &step_xtc, &time_xtc,
                        box_xtc, x_xtc, &prec_xtc);
  if (result_xtc != 0) {
    NOTE("reached the end of XTC file " << fileNameStr);
    endXTC = 1;
  }
  for (int i = 0; i < natom(); ++i) {
    for (int dim = 0; dim < dimen_; ++dim) {
      x_[dimen_*i+dim] = x_xtc[i][dim];
    }
  }

  free(x_xtc);
  delete [] fn_xtc;
  return endXTC;
}
#endif  // XDRFILE_H_

#ifdef XDRFILE_H_
void Space::writeXTC(XDRFILE* trjFileXDR) {
  int natoms_xtc = natom();
  matrix box_xtc;
  box_xtc[0][0] = l_[0];
  box_xtc[0][1] = 0;
  box_xtc[0][2] = 0;
  box_xtc[1][0] = l_[1];
  box_xtc[1][1] = 0;
  box_xtc[1][2] = 0;
  box_xtc[2][0] = l_[2];
  box_xtc[2][1] = 0;
  box_xtc[2][2] = 0;
  // rvec *x_xtc;
  // x_xtc = (rvec *)malloc(sizeof(rvec)*(natom()+1));
  rvec *x_xtc;
  x_xtc = reinterpret_cast<rvec *>(calloc(natoms_xtc, sizeof(x_xtc[0])));
  float prec_xtc = 1000.0;
  for (int i = 0; i < natoms_xtc; ++i) {
    for (int dim = 0; dim < dimen_; ++dim) {
      x_xtc[i][dim] = x(i, dim);
//      cout << "i " << i << " dim " << dim << " x " << x_xtc[i][dim] << endl;
    }
  }
  if (write_xtc(trjFileXDR, natoms_xtc, 0, 0, box_xtc, x_xtc, prec_xtc) != 0) {
    ASSERT(0, "error writing xtc file");
  }
  free(x_xtc);
}
#endif  // XDRFILE_H_

void Space::wrapMol() {
  for (int i = 0; i < nMol(); ++i) {
    wrap(imol2mpart(i));
  }
}

void Space::setAtomAsCOM(const int atom, const vector<int> mpart) {
  vector<double> xnew(3);
  for (unsigned int i = 0; i < mpart.size(); ++i) {
    if (mpart[i] != atom) {
      for (int dim = 0; dim < dimen_; ++dim) {
        xnew[dim] += x(mpart[i], dim) / static_cast<double>(mpart.size()-1);
      }
    }
  }
  for (int dim = 0; dim < dimen_; ++dim) {
    xset(xnew[dim], atom, dim);
  }
}

void Space::setAtomInSphere(const int iAtom, const int jAtom, const double r) {
  vector<double> xnew = ranShell(r, r, dimen_);
  for (int dim = 0; dim < dimen_; ++dim) {
    xset(x(jAtom, dim)+xnew[dim], iAtom, dim);
  }
  if (cellType_ != 0) updateCellofiMol(mol_[iAtom]);
}

void Space::setAtomInCircle(
  const int iAtom,
  const int jAtom,
  const int kAtom,
  const double r,
  const double theta) {
  // to begin, define xnew with origin on the previous atom,
  // and vector pointing in the direction of the bond between previous
  // two atoms, termed AB
  vector<double> xnew(dimen_), ab(dimen_);
  double rsq = 0.;
  for (int dim = 0; dim < dimen_; ++dim) {
    xnew[dim] = x(jAtom, dim) - x(kAtom, dim);
    rsq += pow(xnew[dim], 2.);
  }
  for (int dim = 0; dim < dimen_; ++dim) {
    ab[dim] = xnew[dim]/sqrt(rsq);
    xnew[dim] = ab[dim] * r;
  }

  // rotate xnew by angle (PI-theta) along arbitrary axis which
  // is perpendicular to AB
  vector<double> ortho = orthogonalVec(ab);
  xnew = rotateVecByAxisAngle(xnew, ortho, (PI-theta));

  // rotate xnew by random angle [0, 2PI] along AB
  const double ranAngle = 2*PI*uniformRanNum();
  xnew = rotateVecByAxisAngle(xnew, ab, ranAngle);

  // finally, translate back to original frame of reference
  for (int dim = 0; dim < dimen_; ++dim) {
    xset(x(jAtom, dim)+xnew[dim], iAtom, dim);
  }

  if (cellType_ != 0) updateCellofiMol(mol_[iAtom]);
}

void Space::modBondAngle(const int iAtom, const int jAtom, const int kAtom,
  const double theta) {
  // find the vector orthogonal, ortho, to the vector IJ and KJ
  //  use frame of reference where jAtom is on the origin
  vector<double> xnew(dimen_), ij(dimen_), kj(dimen_);
  for (int dim = 0; dim < dimen_; ++dim) {
    ij[dim] = x(iAtom, dim) - x(jAtom, dim);
    xnew[dim] = kj[dim] = x(kAtom, dim) - x(jAtom, dim);
  }
  normalizeVec(&ij);
  normalizeVec(&kj);
  vector<double> ortho = crossProd(ij, kj);
  if (fabs(vecDotProd(ortho, ortho)) > 10*doubleTolerance) {
    normalizeVec(&ortho);
    for (int dim = 0; dim < dimen_; ++dim) {
      ortho[dim] *= -1;
    }
  } else {
    ortho = orthogonalVec(kj);
  }

  // place iAtom on kAtom, then rotate the IJ vector along OV axis
  // by angle theta
  xnew = rotateVecByAxisAngle(xnew, ortho, theta);

  // translate back to original frame of reference
  for (int dim = 0; dim < dimen_; ++dim) {
    xset(x(jAtom, dim)+xnew[dim], iAtom, dim);
  }

  // update cell list
  if (cellType_ != 0) updateCellofiMol(mol_[iAtom]);

  // update qMol_ and xMolRef_ (but not xMol_)
  qMolInit(mol_[iAtom]);
}

void Space::modBondAngle(const int angleType, const double theta,
  const char* molType) {
  // update angle parameters
  string molTypeStr(molType);
  shared_ptr<Space> s = findAddMolInList(molTypeStr);
  s->qMolInit();
  s->modAngleParams(angleType, 1, theta);

  // obtain atoms describing angle, <ijk
  vector<int> iAtom, jAtom, kAtom;
  const vector<vector<int> > angleList = s->angleList();
  for (unsigned int iAngle = 0; iAngle < angleList.size(); ++iAngle) {
    if (angleList[iAngle][0] == angleType) {
      iAtom.push_back(angleList[iAngle][1]);
      jAtom.push_back(angleList[iAngle][2]);
      kAtom.push_back(angleList[iAngle][3]);

      // update geometry of the molecule stored in addMolList
      s->modBondAngle(iAtom.back(), jAtom.back(), kAtom.back(), theta);
    }
  }
  s->qMolInit();

  // loop through all molecules of molType, call modBondAngle(i,j,k,theta)
  //  randomly swap i and k
  for (int iMol = 0; iMol < nMol(); ++iMol) {
    if (moltype_[iMol].compare(molType) == 0) {
      for (unsigned int iAngle = 0; iAngle < iAtom.size(); ++iAngle) {
        const int firstAtom = mol2part_[iMol];
        const int jAtomTmp = firstAtom + jAtom[iAngle];
        int iAtomTmp = firstAtom + iAtom[iAngle],
            kAtomTmp = firstAtom + kAtom[iAngle];
        if (uniformRanNum() < 0.5) {
          iAtomTmp = firstAtom + kAtom[iAngle];
          kAtomTmp = firstAtom + iAtom[iAngle];
        }
        modBondAngle(iAtomTmp, jAtomTmp, kAtomTmp, theta);
      }
    }
  }
}

/* For coordinates of atom 3, x,y,z, letting atom 4 be origin, solve
 * for three eq and 3 unknowns
 * x^2+y^2+z^2 = L^2
 * x*x1+y*y1+z*z1 = cost143   **x1,x2,etc, are unit normal vectors
 * along 14 or 24 bonds
 * x*x2+y*y2+z*z2 = cost243
 * solve for (two values of) x by substitution and quadratic eq. pick
 * one solution randomly
 * if y1 != 0, H = z2 - z1*y2/y1
 * if H != 0, A = (x1y2/y1 - x2)/H, B = (cost243-cost143*y2/y1)/H
 * z(x) = A*x+B, y(x) = C*x+D
 * C = -x1/y1 - Az1/y1, D = cost143/y1 - Bz1/y1
 * [1+C^2+A^2] x^2 + [2CD+2AB] x + [D^2+B^2-L^2] = 0
 *
 * this solution is plagued by numerical stability, if y1 ~ 0 or |H|<1e-8
 *   modified to do alternative solves for more stable H  */
void Space::setAtomInBranch(const int a1, const int a2, const int a3,
  const int a4, const double t143, const double t243, const double L) {
  ASSERT(dimen_ == 3, "setAtomInBranch must have dimen(" << dimen_ << ") == 3");

  // let atom 4 be origin, define vectors r1, r2 as unit normals
  // along 14 and 24 bonds
  double x1, y1, z1, x2, y2, z2, r;
  x1 = x(a1, 0) - x(a4, 0);
  y1 = x(a1, 1) - x(a4, 1);
  z1 = x(a1, 2) - x(a4, 2);
  x2 = x(a2, 0) - x(a4, 0);
  y2 = x(a2, 1) - x(a4, 1);
  z2 = x(a2, 2) - x(a4, 2);
  r = sqrt(x1*x1+y1*y1+z1*z1);
  x1 /= r; y1 /= r; z1 /= r;
  r = sqrt(x2*x2+y2*y2+z2*z2);
  x2 /= r; y2 /= r; z2 /= r;
  const double c143 = cos(t143), c243 = cos(t243);
  double x3, y3, z3;
  if ( fabs(y1) > fabs(x1) ) {
//  if ( (fabs(x1) < doubleTolerance) || (fabs(Cyz) > fabs(Cxz)) ) {
    // cout << "test1 " << fabs(x1) << " t2 " << fabs(Cyz) << " > "
    //      << fabs(Cxz) << endl;
    solveBranch_(x1, y1, z1, x2, y2, z2, &x3, &y3, &z3, c143, c243);
  } else {
  // } else if ( (fabs(y1) < doubleTolerance) || (fabs(Cyz) < fabs(Cxz)) ) {
    // cout << "test2 " << fabs(y1) << " t2 " << fabs(Cyz) << " > "
    //      << fabs(Cxz) << endl;
    solveBranch_(y1, x1, z1, y2, x2, z2, &y3, &x3, &z3, c143, c243);
  }
  xset(L*x3+x(a4, 0), a3, 0);
  xset(L*y3+x(a4, 1), a3, 1);
  xset(L*z3+x(a4, 2), a3, 2);
  if (cellType_ != 0) updateCellofiMol(mol_[a1]);
}

vector<double> Space::bondParams(const int iAtom, const int jAtom) {
  double k, l0;
  ASSERT(bondList_.size() > 0, "no bonds when searching for bondParams");

  int i = 0;
  while ( ( ( (bondList_[i][1] != iAtom) || (bondList_[i][2] != jAtom) ) &&
            ( (bondList_[i][1] != jAtom) || (bondList_[i][2] != iAtom) ) ) &&
          (i < static_cast<int>(bondList_.size()) ) ) {
    ++i;
  }
  if (i != static_cast<int>(bondList_.size())) {
    const int btype = bondList_[i][0];
    k  = bondParam_[btype][0];
    l0 = bondParam_[btype][1];
  }
  vector<double> params;
  params.push_back(k);
  params.push_back(l0);
  return params;
}

vector<double> Space::angleParams(const int iAtom, const int jAtom,
  const int kAtom) {
  ASSERT(angleList_.size() > 0, "no angles when searching for angleParams");
  int i = 0;
  while ( ( ( (angleList_[i][1] != iAtom) ||
              (angleList_[i][2] != jAtom) ||
              (angleList_[i][3] != kAtom) ) &&
            ( (angleList_[i][1] != iAtom) ||
              (angleList_[i][2] != kAtom) ||
              (angleList_[i][3] != jAtom) ) &&
            ( (angleList_[i][1] != jAtom) ||
              (angleList_[i][2] != iAtom) ||
              (angleList_[i][3] != kAtom) ) &&
            ( (angleList_[i][1] != jAtom) ||
              (angleList_[i][2] != kAtom) ||
              (angleList_[i][3] != iAtom) ) &&
            ( (angleList_[i][1] != kAtom) ||
              (angleList_[i][2] != iAtom) ||
              (angleList_[i][3] != jAtom) ) &&
            ( (angleList_[i][1] != kAtom) ||
              (angleList_[i][2] != jAtom) ||
              (angleList_[i][3] != iAtom) ) ) &&
          (i < static_cast<int>(angleList_.size()) ) ) {
    ++i;
  }
  vector<double> params;
  if (i != static_cast<int>(angleList_.size())) {
    const int atype = angleList_[i][0];
    params.push_back(angleParam_[atype][0]);
    params.push_back(angleParam_[atype][1]);
  }
  return params;
}

vector<vector<int> > Space::listBonds(const int iAtom) {
  vector<vector<int> > list;
  const int iMol = mol_[iAtom];
  const int ia = iAtom - mol2part_[iMol];
  shared_ptr<Space> s = findAddMolInList(moltype_[iMol]);
  vector<vector<int> > bList = s->bondList();
  for (unsigned int i = 0; i < bList.size(); ++i) {
    if ( (bList[i][1] == ia) || (bList[i][2] == ia) ) {
      list.push_back(bList[i]);
    }
  }
  return list;
}

vector<vector<int> > Space::listAngles(const int iAtom, const int jAtom) {
  vector<vector<int> > list;
  const int iMol = mol_[iAtom];
  const int ia = iAtom - mol2part_[iMol];
  const int ja = jAtom - mol2part_[iMol];
  shared_ptr<Space> s = findAddMolInList(moltype_[iMol]);
  vector<vector<int> > aList = s->angleList();
  for (unsigned int i = 0; i < aList.size(); ++i) {
    if ( ( (aList[i][1] == ia) && (aList[i][2] == ja) ) ||
         ( (aList[i][1] == ia) && (aList[i][3] == ja) ) ||
         ( (aList[i][2] == ia) && (aList[i][1] == ja) ) ||
         ( (aList[i][2] == ia) && (aList[i][3] == ja) ) ||
         ( (aList[i][3] == ia) && (aList[i][1] == ja) ) ||
         ( (aList[i][3] == ia) && (aList[i][2] == ja) ) ) {
      list.push_back(aList[i]);
    }
  }
  return list;
}

void Space::solveBranch_(const double x1, const double y1, const double z1,
  const double x2, const double y2, const double z2, double *x3, double *y3,
  double *z3, const double c143, const double c243) {
  ASSERT(y1 != 0, "y1==0");
  const long double H = z2 - y2*z1/y1;
  ASSERT(H != 0, "H==0");
  const long double A = (x1*y2/y1 - x2)/H,
               B = (c243 - c143*y2/y1)/H,
               C = -x1/y1 - A*z1/y1,
               D = c143/y1 - B*z1/y1,
               a = (1+A*A+C*C),
               b = 2*(A*B+C*D),
               c = (B*B+D*D-1);
  long double ans1, ans2;
  const long double discrim = quadraticEqReal(a, b, c, ans1, ans2);
  if (discrim < 0) {
    if ( (sqrt(fabs(discrim))/2/fabs(a) < 1000*sqrt(doubleTolerance)) ||
         fabs(discrim) < 10000*doubleTolerance) {
      // within double preicison, the discriminant is zero
      ans1 = ans2 = -b/2/a;
    } else {
      std::streamsize ss = cout.precision();
      cout << std::setprecision(std::numeric_limits<long double>::digits10+2)
           << "c143 " << c143 << " c243 " << c243 << endl;
      cout << "x1 " << x1 << " " << y1 << " " << z1 << endl;
      cout << "x2 " << x2 << " " << y2 << " " << z2 << endl;
      cout << "A " << A << " H " << H << " B " << B << " C " << C
           << " D " << D << endl;
      cout << "ans1 " << ans1 << " ans2 " << ans2 << endl;
      cout << "discrim " << discrim << endl;
      cout << "a " << a << " b " << b << " c " << c << endl;
      cout << std::setprecision(ss);
      cout << "tol " << sqrt(doubleTolerance) << " rel "
           << sqrt(fabs(discrim))/2/fabs(a) << endl;
      cout << "tol " << 10*doubleTolerance << " rel " << fabs(discrim)
           << endl;
      ASSERT(0, "imaginary branch");
    }
  }
  if (uniformRanNum() < 0.5) {
    ans1 = ans2;
  }

  // (x,y,z) is unit vector pointing in direction of l3 bond
  *x3 = ans1;
  *y3 = C*ans1+D;
  *z3 = A*ans1+B;
}

void Space::nRadialHist(Histogram *nhistPtr) {
  Histogram& nhist = *nhistPtr;
  nhist.count();
  const int itype = nhist.iType(), jtype = nhist.jType();
  double lx = l_[0], ly = 0, lz = 0, halflx = lx/2.,
    halfly = 0, halflz = 0, dx, dy, dz, xi, yi, zi = 0, r2 = 0;
  const double rCutSq = pow(minl()/2., 2);
  if (dimen_ >= 2) {
    ly = l_[1], halfly = ly/2.;
  }
  if (dimen_ >= 3) {
    lz = l_[2], halflz = lz/2.;
  }
  for (int ipart = 0; ipart < natom(); ++ipart) {
    if (type_[ipart] == itype) {
      if (dimen_ == 2) {
        xi = x_[dimen_*ipart]; yi = x_[dimen_*ipart+1];
      } else {
        xi = x_[dimen_*ipart]; yi = x_[dimen_*ipart+1]; zi = x_[dimen_*ipart+2];
      }
      for (int jpart = 0; jpart < natom(); ++jpart) {
        if (type_[jpart] == jtype) {
          if (mol_[ipart] != mol_[jpart]) {
            if (dimen_ == 2) {
              dx = xi - x_[dimen_*jpart];
              dy = yi - x_[dimen_*jpart+1];
              if (dx >  halflx) dx -= lx;
              if (dx < -halflx) dx += lx;
              if (dy >  halfly) dy -= ly;
              if (dy < -halfly) dy += ly;
              r2 = dx*dx + dy*dy;
            } else if (dimen_ == 3) {
              dx = xi - x_[dimen_*jpart];
              dy = yi - x_[dimen_*jpart+1];
              dz = zi - x_[dimen_*jpart+2];
              if (dx >  halflx) dx -= lx;
              if (dx < -halflx) dx += lx;
              if (dy >  halfly) dy -= ly;
              if (dy < -halfly) dy += ly;
              if (dz >  halflz) dz -= lz;
              if (dz < -halflz) dz += lz;
              r2 = dx*dx + dy*dy + dz*dz;
            }
            if (r2 <= rCutSq) nhist.accumulate(sqrt(r2));
          }
        }
      }
    }
  }
}

void Space::printRadial(const Histogram &nhist, const char* fileName) {
  std::ofstream file(fileName);
  file << "# iType " << nhist.iType() << " jType " << nhist.jType() << endl;
  file << "# r g nhist" << endl;
  file << nhist.bin2m(0) - nhist.binWidth() << " 0 0" << endl;
  for (int i = 0; i < nhist.size(); ++i) {
    const double r = nhist.bin2m(i), rmin = r - 0.5*nhist.binWidth(),
                 rmax = r + 0.5*nhist.binWidth();
    double fac = -1;
    if (dimen() == 3) {
      fac = 4./3.;
    } else if (dimen() == 2) {
      fac = 1;
    } else {
      ASSERT(0, "this dimensionality(" << dimen() << ") is not implemented"
             << "for printing a radial disribution function");
    }
    const double nideal = fac*PI*(nMol()/vol())*
        (pow(rmax, dimen())-pow(rmin, dimen()));
    file << r << " " << static_cast<double>(nhist.hist()[i]) /
      static_cast<double>(nhist.nCount()) /
      static_cast<double>(nMol())/nideal << " " << nhist.hist()[i] << endl;
  }
}

void Space::pivotMol(const int iMol, const vector<double> r) {
  for (int iAtom = mol2part_[iMol]; iAtom < mol2part_[iMol+1]; ++iAtom) {
    for (int dim = 0; dim < dimen_; ++dim) {
      x_[dimen_*iAtom+dim] = 2*r[dim] - x_[dimen_*iAtom+dim];
    }
  }
}

vector<double> Space::randPosition() {
  vector<double> x(dimen_);
  for (int dim = 0; dim < dimen_; ++dim) {
    ASSERT(l_[dim] != 0, "the domain must be set for randPosition, l_[dim="
      << dim << "]=" << l_[dim]);
    x[dim] = l_[dim]*(uniformRanNum() - 0.5);
  }
  return x;
}

vector<double> Space::randPosition(const double iMol, const double maxDisp) {
  vector<double> xtmp(dimen_);
  for (int dim = 0; dim < dimen_; ++dim) {
    ASSERT(l_[dim] != 0, "the domain must be set for randPosition, l_[dim="
      << dim << "]=" << l_[dim]);
    if (maxDisp == -1) {
      xtmp[dim] = l_[dim]*(uniformRanNum() - 0.5);
    } else {
      xtmp[dim] = x(mol2part_[iMol], dim) + maxDisp*(2.*uniformRanNum() - 1.);
    }
  }
  rwrap(&xtmp);
  return xtmp;
}

vector<double> Space::scatterIntensity(const double qMin, const double qMax,
  const double dq) {
  ASSERT(dimen_ == 3,
    "dimen(" << dimen_ << ") must equal 3 for scatterIntensity");
  ASSERT(minl() != 0, "scatterIntensity assumes PBC");

  const int nq = (qMax-qMin)/dq + 1;
  vector<double> intensity(nq);
  double r2, xi, yi, zi, dx, dy, dz;
  const double lx = l_[0], ly = l_[1], lz = l_[2], halflx = lx/2.,
    halfly = ly/2., halflz = lz/2.;
  for (int iAtom = 0; iAtom < natom() - 1; ++iAtom) {
    xi = x(iAtom, 0);
    yi = x(iAtom, 1);
    zi = x(iAtom, 2);
    for (int jAtom = iAtom + 1; jAtom < natom(); ++jAtom) {
      // separation distance with periodic boundary conditions
      dx = xi - x_[dimen_*jAtom];
      dy = yi - x_[dimen_*jAtom+1];
      dz = zi - x_[dimen_*jAtom+2];
      if (dx >  halflx) dx -= lx;
      if (dx < -halflx) dx += lx;
      if (dy >  halfly) dy -= ly;
      if (dy < -halfly) dy += ly;
      if (dz >  halflz) dz -= lz;
      if (dz < -halflz) dz += lz;
      r2 = dx*dx + dy*dy + dz*dz;
      const double r = sqrt(r2);

      for (int iq = 0; iq < nq; ++iq) {
        const double qr = (qMin + iq*dq)*r;
        intensity[iq] += sin(qr)/qr;
      }
    }
  }
  return intensity;
}

void Space::scaleDomain(const double factor, const int dim) {
  ASSERT((dim <= dimen_) && (dim >= 0),
   "dim(" << dim << ") in scaleDomain is outside of range for dimen("
   << dimen_ << ")");
  ASSERT(factor > 0,
    "factor(" << factor << ") in scaleDomain cannot be negative");

  // scale the box subject to bounds
  double factorActual = factor;
  if (maxlFlag_ != 0) {
    // check that the box isn't scaled beyond limits
    const double lNew = l_[dim]*factor;
    if (lNew > maxl_[dim]) {
      factorActual = maxl_[dim]/l_[dim];
    }
  }
  lset(l_[dim]*factorActual, dim);

  // loop through each molecule, and scale based on the position
  //  of the first site in molecule
  //  if molecules are one site, then reduces to scaling the position
  for (int iMol = 0; iMol < nMol(); ++iMol) {
    const double dx = (factorActual - 1.)*x(mol2part_[iMol], dim);
    for (int iAtom = mol2part_[iMol]; iAtom < mol2part_[iMol+1]; ++iAtom) {
      x_[dimen_*iAtom + dim] += dx;
    }
  }

  if (cellType() > 0) updateCells();
}

void Space::avb(const int iAtom, const int jAtom, const double rAbove,
                const double rBelow, const char* region) {
  string regionStr(region);
  vector<double> xnew;
  if (regionStr.compare("in") == 0) {
    // generate xnew as random position in shell,
    // then shift by position of jAtom
    xnew = ranShell(rAbove, rBelow, dimen_);
    for (int dim = 0; dim < dimen_; ++dim) {
      xnew[dim] += x(jAtom, dim);
    }

  } else if (regionStr.compare("out") == 0) {
    // pick random position in box, reject if in bonded region of jAtom
    int term = 0, n = 0, nMax = 1e6;
    vector<double> xj;
    for (int dim = 0; dim < dimen_; ++dim) xj.push_back(x_[dimen_*jAtom+dim]);
    while ( (term == 0) && (n < nMax) ) {
      xnew = randPosition();
      const double r2 = rsq(xnew, xj);
      if ( (r2 > rAbove*rAbove) || (r2 < rBelow*rBelow) ) term = 1;
      ++n;
    }
    ASSERT(n < nMax,
      "reached limiting number of attempts, nMax(" << nMax << ") in avb out");

  } else {
    ASSERT(0, "unrecognized region(" << region << ") in avb");
  }

  xset(iAtom, xnew);
}

double Space::Q6(const double rCut
  ) {
  xMolGen();
  vector<double> rij(3);
  vector<std::complex<double> > sphH(13);
  int nBond = 0;
  for (int iMol = 0; iMol < nMol()-1; ++iMol) {
    const vector<double> xi = xMol()[iMol][0];
    for (int jMol = iMol+1; jMol < nMol(); ++jMol) {
      const vector<double> xj = xMol()[jMol][0];
      for (int dim = 0; dim < dimen(); ++dim) {
        rij[dim] = xi[dim] - xj[dim];
      }
      const vector<double> dx = pbc(rij);
      for (int dim = 0; dim < dimen(); ++dim) {
        rij[dim] += dx[dim];
      }
      const double r2 = vecDotProd(rij, rij);
      if (r2 < rCut*rCut) {
        const vector<std::complex<double> > sphHTmp = cart2sphereHarm6(rij);
        for (unsigned int i = 0; i < sphH.size(); ++i) {
          sphH[i] += sphHTmp[i];
        }
        ++nBond;
      }
    }
  }
  for (unsigned int i = 0; i < sphH.size(); ++i) {
    sphH[i] /= static_cast<double>(nBond);
  }
  return sqrt(4*PI/13.*complexVec2norm(sphH));
}

void Space::setXYTilt(const double xyTilt) {
  ASSERT(xyTilt <= l_[0], "the xyTilt(" << xyTilt << ") cannot be"
    << "larger than the box(" << l_[0] << ")");
  xyTilt_ = xyTilt;
}

void Space::setXZTilt(const double xzTilt) {
  ASSERT(xzTilt <= l_[0], "the xzTilt(" << xzTilt << ") cannot be"
    << "larger than the box(" << l_[0] << ")");
  xzTilt_ = xzTilt;
}

void Space::setYZTilt(const double yzTilt) {
  ASSERT(yzTilt <= l_[1], "the yzTilt(" << yzTilt << ") cannot be"
    << "larger than the box(" << l_[1] << ")");
  yzTilt_ = yzTilt;
}

void Space::modXYTilt(const double deltaXYTilt) {
  floppyBox_ = 1;
  const double xyTiltOld = xyTilt_;
  xyTilt_ += deltaXYTilt;

  // limit xyTilt to some percentage of the box
  const double maxPercBox = 0.25;
  if (xyTilt_ >  maxPercBox*l_[0]) xyTilt_ =  maxPercBox*l_[0];
  if (xyTilt_ < -maxPercBox*l_[0]) xyTilt_ = -maxPercBox*l_[0];

  // transform the particles, based on the position of the first site
  // of each 'molecule'
  for (int iMol = 0; iMol < nMol(); ++iMol) {
    const double dx = x(mol2part_[iMol], 1) * (xyTilt_ - xyTiltOld) / l_[1];

    for (int iAtom = mol2part_[iMol]; iAtom < mol2part_[iMol+1]; ++iAtom) {
      xset(x(iAtom, 0) + dx, iAtom, 0);
    }
  }
}

void Space::modXZTilt(const double deltaXZTilt) {
  floppyBox_ = 1;
  const double xzTiltOld = xzTilt_;
  xzTilt_ += deltaXZTilt;

  // limit xzTilt to some percentage of the box
  const double maxPercBox = 0.25;
  if (xzTilt_ >  maxPercBox*l_[0]) xzTilt_ =  maxPercBox*l_[0];
  if (xzTilt_ < -maxPercBox*l_[0]) xzTilt_ = -maxPercBox*l_[0];

  // transform the particles, based on the position of the first site
  // of each 'molecule'
  for (int iMol = 0; iMol < nMol(); ++iMol) {
    const double dx = x(mol2part_[iMol], 2) * (xzTilt_ - xzTiltOld) / l_[2];

    for (int iAtom = mol2part_[iMol]; iAtom < mol2part_[iMol+1]; ++iAtom) {
      xset(x(iAtom, 0) + dx, iAtom, 0);
    }
  }
}

void Space::modYZTilt(const double deltaYZTilt) {
  floppyBox_ = 1;
  const double yzTiltOld = yzTilt_;
  yzTilt_ += deltaYZTilt;

  // limit yzTilt to some percentage of the box
  const double maxPercBox = 0.25;
  if (yzTilt_ >  maxPercBox*l_[1]) yzTilt_ =  maxPercBox*l_[1];
  if (yzTilt_ < -maxPercBox*l_[1]) yzTilt_ = -maxPercBox*l_[1];

  // transform the particles, based on the position of the first site
  // of each 'molecule'
  for (int iMol = 0; iMol < nMol(); ++iMol) {
    const double dx = x(mol2part_[iMol], 2) * (yzTilt_ - yzTiltOld) / l_[2];

    for (int iAtom = mol2part_[iMol]; iAtom < mol2part_[iMol+1]; ++iAtom) {
      xset(x(iAtom, 1) + dx, iAtom, 1);
    }
  }
}

double Space::minBondLength() {
  double min = 1e50;
  xMolGen();

  // loop through all existing molecules
  for (unsigned int iMol = 0; iMol < xMol_.size(); ++iMol) {
    for (unsigned int iPart = 0; iPart < xMol_[iMol].size(); ++iPart) {
      for (unsigned int jPart = 0; jPart < xMol_[iMol].size(); ++jPart) {
        if (iPart != jPart) {
          const double r2 = rsq(xMol_[iMol][iPart], xMol_[iMol][jPart]);
          if (r2 < min) min = r2;
        }
      }
    }
  }

  // loop through potential molecules to add
  for (unsigned int iAdd = 0; iAdd < addMolList_.size(); ++iAdd) {
    const double r2 = pow(addMolList_[iAdd]->minBondLength(), 2.);
    if (r2 < min) min = r2;
  }
  return sqrt(min);
}

/**
 * return quaternion vector of molecule
 */
vector<double> Space::qMol(const int iMol) const {
  if (eulerFlag_ == 0) {
    vector<double> q(qdim_);
    for (int dim = 0; dim < qdim_; ++dim) {
      q[dim] = qMol_[qdim_*iMol + dim];
    }
    return q;
  } else {
    vector<double> q(dimen_);
    for (int dim = 0; dim < dimen_; ++dim) {
      q[dim] = qMol_[qdim_*iMol + dim];
    }
    return q;
  }
}

// **
// * update the euler angles of iMol according to position relative to ref
// */
// void Space::pos2euler(const int iMol) {
//  vector<vector<double> > xmol;
//  for (int iatom = mol2part_[iMol]; iatom < mol2part_[iMol+1]; ++iatom) {
//    vector<double> xa(dimen_);
//    for (int dim = 0; dim < dimen_; ++dim) {
//      xa[dim] = x(iatom, dim);
//    }
//    xmol.push_back(xa);
//  }
//  cout << "q " << qMol_[0] << " " << qMol_[1] << endl;
//  vec2str(xMolRef_[iMol]);
//  cout << "xMolRef " << xMolRef_[iMol][0][0] << endl;
//  vector<vector<double> > xrefinv = inv3by3(xMolRef_[iMol]);
//  cout << "xrefinv " << xrefinv[0][0] << endl;
//  vector<vector<double> > rotMat = matMul(xmol, xrefinv);
//  cout << "rots " << rotMat.size() << " " << rotMat[0].size() << endl;
//  cout << rotMat[0][0] << " " << rotMat[0][1] << " " << rotMat[0][2] << endl;
//  cout << rotMat[1][0] << " " << rotMat[1][1] << " " << rotMat[1][2] << endl;
//  cout << rotMat[2][0] << " " << rotMat[2][1] << " " << rotMat[2][2] << endl;
//  vector<vector<double> > euler = RotMat2Euler(rotMat);
//  exit(0);
// }
//
//

int Space::randMolofType(const int iType) {
  int iMol = -1, iter = 0, max = 10*nMol();
  bool term = false;
  ASSERT((iType >= 0) && (iType < static_cast<int>(addMolList().size())),
    "Invalid iType supplied");
  while (!term && (iter < max)) {
    iMol = uniformRanNum(0, nMol()-1);
    if (addMolListType_[iType].compare(moltype_[iMol]) == 0) {
      term = true;
    }
    ++iter;
  }
  return iMol;
}

void Space::swapPositions(const int iMol, const int jMol) {
  ASSERT(nMol() == natom(), "swap implemented for monoatomics");
  for (int dim = 0; dim < dimen_; ++dim) {
    const double dx = x_[dimen_*mol2part_[iMol]+dim]
                    - x_[dimen_*mol2part_[jMol]+dim];
    for (int iatom = mol2part_[iMol]; iatom < mol2part_[iMol+1]; ++iatom) {
      x_[dimen_*iatom+dim] -= dx;
    }
    for (int iatom = mol2part_[jMol]; iatom < mol2part_[jMol+1]; ++iatom) {
      x_[dimen_*iatom+dim] += dx;
    }
  }
  if (cellType_ > 0) updateCellofiMol(iMol);
  if (cellType_ > 0) updateCellofiMol(jMol);
}

void Space::printxyzvmd(const char* fileName, const int initFlag) {
  // write vmd script to visualize in fileName appended with .vmd
  if ( (initFlag == 1) || (initFlag == 2) ) {
    std::stringstream vmdfnamess;
    vmdfnamess << fileName << ".vmd";
    std::ofstream vmdf(vmdfnamess.str().c_str());
    double radius;
    vmdf << "display projection Orthographic" << endl
      << "color Display Background white" << endl
      << "axes location Off" << endl;
    if (initFlag == 2) {
      radius = 10;
      vmdf << "mol load xyz " << trim("/", fileName) << ".xyz" << endl
           << "animate delete beg 0 end 0" << endl
           << "mol addfile " << trim("/", fileName) << ".xtc" << endl;
    } else {
      radius = 1;
      vmdf << "topo readvarxyz " << trim("/", fileName) << ".xyz" << endl;
    }
    vmdf << "mol modstyle 0 0 VDW 1.0000000 120.000000" << endl
      << "set sel [atomselect top all]" << endl
      << "$sel set radius " << 0.5*radius << endl
      << "$sel set mass 1" << endl;
  }
}

vector<double> Space::ipart2euler(const int ipart) {
  vector<double> xvec(dimen_), xref(dimen_, 0.);
  xref[dimen_-1] = 1.;
  for (int dim = 0; dim < dimen_; ++dim) {
    xvec[dim] = x_[dimen_*(ipart+1)+dim] - x_[dimen_*ipart+dim];
  }
  vector<vector<double> > a1 = outerProd(xvec, xref),
    a2 = outerProd(xref, xref),
    a3 = transpose(a2),
    rot = matMul(a1, a3),
    eulertmp = RotMat2Euler(rot);
  return eulertmp[0];
}


#ifdef JSON_
/**
 * Initialize with JSON data file
 * Reads number of atoms, molecules, atom types, masses
 * Resizes appropriate arrays
 */
void Space::initJSONData(
  const std::string fileName,  //!< LAMMPS Data file name
  const int nTypesExist  //!< number of particle types that already exists
  ) {
  // open a JSON file
  std::ifstream file(fileName.c_str());
  ASSERT(file.good(), "cannot find json DATA file " << fileName);
  nlohmann::json j;
  file >> j;

  ASSERT(j["dimen"] == dimen_, "dimensions of space class(" << dimen_ <<
    ") do not match data file(" << j["dimen"] << ") " << fileName);

  // atom properties
  mol2part_.clear();
  listAtoms_.clear();
  nType_.clear();
  int imolprev = 0;
  for (unsigned int i = 0; i < j["atoms"].size(); ++i) {
    for (int dim = 0; dim < dimen_; ++dim) {
      stringstream ss;
      ss << "dim" << dim;
      const double x = j["atoms"][i][ss.str().c_str()];
      x_.push_back(x);
    }
    const int imol = j["atoms"][i]["mol"];
    mol_.push_back(imol - 1);
    const int itype = j["atoms"][i]["type"];
    type_.push_back(itype - 1 + nTypesExist);
    if (static_cast<int>(nType_.size()) <= itype - 1) {
      nType_.resize(itype);
    }
    ++nType_[itype - 1];
    listAtoms_.push_back(i);

    // new molecule?
    if ( (i == 0) || (imol != imolprev) ) {
      imolprev = imol;

      // add new molecule
      mol2part_.push_back(i);
      moltype_.push_back(fileName);
    }
  }
  mol2part_.push_back(natom());
  xMolGen();
  if (natom() != 1) {
    qMolInit();
  }

  // bond properties
  string typestr("bonds");
  for (int index = 0; index < static_cast<int>(j[typestr].size()); ++index) {
    vector<int> b;
    const int type = j[typestr][index]["type"]; b.push_back(type-1);
    const int a0   = j[typestr][index]["a0"];   b.push_back(a0-1);
    const int a1   = j[typestr][index]["a1"];   b.push_back(a1-1);
    bondList_.push_back(b);
  }
  typestr.assign("bondCoeffs");
  for (int ib = 0; ib < static_cast<int>(j[typestr].size()); ++ib) {
    vector<double> b;
    const double k = j[typestr][ib]["k"]; b.push_back(k);
    const double L = j[typestr][ib]["L"]; b.push_back(L);
    bondParam_.push_back(b);
  }

  // angle properties
  typestr.assign("angles");
  for (int index = 0; index < static_cast<int>(j[typestr].size()); ++index) {
    vector<int> b;
    const int type = j[typestr][index]["type"]; b.push_back(type-1);
    const int a0   = j[typestr][index]["a0"];   b.push_back(a0-1);
    const int a1   = j[typestr][index]["a1"];   b.push_back(a1-1);
    const int a2   = j[typestr][index]["a2"];   b.push_back(a2-1);
    angleList_.push_back(b);
  }
  typestr.assign("angleCoeffs");
  for (int ib = 0; ib < static_cast<int>(j[typestr].size()); ++ib) {
    vector<double> b;
    const double k = j[typestr][ib]["k"]; b.push_back(k);
    const double theta = j[typestr][ib]["theta"]; b.push_back(theta/180*PI);
    angleParam_.push_back(b);
  }
}
#endif  // JSON_

void Space::initData(const std::string fileName, const int nTypesExist) {
  // use file extension to determine whether to use JSON or LMP data files
  if (trim(".", fileName) == "json") {
    #ifdef JSON_
      initJSONData(fileName, nTypesExist);
    #else
      ASSERT(0, "JSON files are not implemented, but file " << fileName <<
        "as attempted to be read");
    #endif  // JSON_

  } else {
    initLMPData(fileName, nTypesExist);
  }
}

void Space::initIntra(const vector<vector<int> >& map) {
  const int iMolType = nMolTypes() - 1;
  intraMap_.resize(nParticleTypes(), vector<vector<int> >(
    map.size(), vector<int>(map.size())));
  for (int im = 0; im < static_cast<int>(map.size()); ++im) {
    for (int jm = 0; jm < static_cast<int>(map[im].size()); ++jm) {
      intraMap_[iMolType][im][jm] = map[im][jm];
    }
  }
}

void Space::replicate(const int nx, const int ny, const int nz) {
  // catch implementation caveat
  ASSERT(dimen() == 3, "replicate assumes 3D");
  ASSERT( (xyTilt_ == 0) && (xzTilt_ == 0) && (yzTilt_ == 0),
    "tilt is not implemented in replication");

  // store original variables
  int ipartBig = natom(), iMolBig = nMol();
  const int nMolOrig = nMol();
  const vector<double> boxOrig = l_;

  for (int dim = 0; dim < dimen(); ++dim) {
    lset(2*l_[dim], dim);
  }
  for (int ix = 0; ix <= nx; ++ix) {
  for (int iy = 0; iy <= ny; ++iy) {
  for (int iz = 0; iz <= nz; ++iz) {
    if ( (ix != 0) || (iy != 0) || (iz != 0) ){
      for (int iMol = 0; iMol < nMolOrig; ++iMol) {
        vector<int> lshift;
        lshift.push_back(ix);
        lshift.push_back(iy);
        lshift.push_back(iz);
        addMol(moltype_[iMol].c_str());
        for (int ipart = mol2part_[iMol]; ipart < mol2part_[iMol+1]; ++ipart) {
          for (int dim = 0; dim < dimen(); ++dim) {
            // positions
            xset(x(ipart, dim) + lshift[dim]*boxOrig[dim], ipartBig, dim);

            // reference positions
            if (!sphereSymMol_) {
              const int iAtom = ipart - mol2part_[iMol];
              xMolRef_[iMolBig][iAtom][dim] = xMolRef_[iMol][iAtom][dim];
            }
          }
          ++ipartBig;
        }

        // orientations
        if (!sphereSymMol_) {
          for (int qd = 0; qd < qdim_; ++qd) {
            qMolAlt(iMolBig, qd, qMol(iMol, qd));
          }
        }
        ++iMolBig;
      }
    }
  }}}
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_


