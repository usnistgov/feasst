/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./analyze_scatter.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

AnalyzeScatter::AnalyzeScatter(Pair *pair, const argtype &args)
  : Analyze(pair, args) {
  defaultConstruction_();
  argparse_.initArgs(className_, args);
  argparse_.checkAllArgsUsed();
}

AnalyzeScatter::AnalyzeScatter(
  Pair *pair,
  const char* fileName)
  : Analyze(pair, fileName) {
  defaultConstruction_();
  double dgrTmp = fstod("dgrRadialBinDist", fileName);
  int nMacrosTmp = fstoi("nMacros", fileName);
  std::string strtmp = fstos("nMomentsCut", fileName);
  if (!strtmp.empty()) {
    nMomentsCut_ = stoi(strtmp);
  }

  // cout << "nm " << nMacrosTmp << endl;
  initSANS(dgrTmp, nMacrosTmp);

  // cout << " open file and skip header lines" << endl;
  std::ifstream fs(fileName);
  string line;
  const int nLines = numLines(fileName);
  int nSkip = nLines - nMacrosTmp*nbins_;
  for (int i = 0; i < nSkip; ++i) getline(fs, line);

  for (unsigned int iMacro = 0; iMacro < countConf_.size(); ++iMacro) {
    // obtain countConf from commented lines
    stringstream ss;
    ss << "countConf" << iMacro;
    countConf_[iMacro] = fstoll(ss.str().c_str(), fileName);
    for (int bin = 0; bin < nbins_; ++bin) {
      double r;
      fs >> r;
      for (unsigned int iType = 0; iType < histInter2_[iMacro].size();
           ++iType) {
        for (unsigned int jType = iType; jType < histInter2_[iMacro][0].size();
             ++jType) {
          int hval;
          fs >> hval;
          histInter2_[iMacro][iType][jType][bin] = hval;
          histInter2_[iMacro][jType][iType][bin] =
            histInter2_[iMacro][iType][jType][bin];
        }
      }
      for (unsigned int iType = 0; iType < histIntra2_[iMacro].size();
           ++iType) {
        for (unsigned int jType = iType; jType < histIntra2_[iMacro][0].size();
             ++jType) {
          fs >> histIntra2_[iMacro][iType][jType][bin];
          histIntra2_[iMacro][jType][iType][bin] =
            histIntra2_[iMacro][iType][jType][bin];
        }
      }
      for (unsigned int iType = 0; iType < histMoments_.size(); ++iType) {
        for (unsigned int jType = iType; jType < histMoments_[0].size();
             ++jType) {
          for (int iMo = iType; iMo < nMomentsCut_; ++iMo) {
            fs >> histMoments_[iMacro][iType][jType][iMo][bin];
            histMoments_[iMacro][jType][iType][iMo][bin] =
              histMoments_[iMacro][iType][jType][iMo][bin];
          }
        }
      }

      getline(fs, line);
    }
  }
}

void AnalyzeScatter::writeRestart(const char* fileName, const int iMacro) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  file << "# dgrRadialBinDist " << dgr_ << endl;
  file << "# nMacros " << countConf_.size() << endl;
  if (nMomentsCut_ != 0) file << "# nMomentsCut " << nMomentsCut_ << endl;
  for (unsigned int im = 0; im < countConf_.size(); ++im) {
    if ( (iMacro == -1) || (static_cast<int>(im) == iMacro) ) {
      file << "# countConf" << im << " " << countConf_[im] << endl;
    }
  }
  for (unsigned int im = 0; im < countConf_.size(); ++im) {
    if ( (iMacro == -1) || (static_cast<int>(im) == iMacro) ) {
      file << "#r ";
      for (unsigned int iType = 0; iType < histInter2_[im].size(); ++iType) {
        for (unsigned int jType = iType; jType < histInter2_[im][0].size(); ++jType) {
          file << "interi" << iType << "j" << jType << " ";
        }
      }
      for (unsigned int iType = 0; iType < histIntra2_[im].size(); ++iType) {
        for (unsigned int jType = iType; jType < histIntra2_[im][0].size(); ++jType) {
          file << "intrai" << iType << "j" << jType << " ";
        }
      }
      for (unsigned int iType = 0; iType < histMoments_[im].size(); ++iType) {
        for (unsigned int jType = iType; jType < histMoments_[im][0].size(); ++jType) {
          for (int iMo = 0; iMo < nMomentsCut_; ++iMo) {
            file << "hi" << iType << "j" << jType << "m" << iMo << " ";
          }
        }
      }
      file << endl;
      vector<vector<vector<long long> > > &hist = histInter2_[im];
      for (int bin = 0; bin < nbins_; ++bin) {
        const double r = dgr_*(bin + 0.5);
        file << r << " ";
        for (unsigned int iType = 0; iType < hist.size(); ++iType) {
          for (unsigned int jType = iType; jType < hist[0].size(); ++jType) {
            file << histInter2_[im][iType][jType][bin] << " ";
          }
        }
        for (unsigned int iType = 0; iType < hist.size(); ++iType) {
          for (unsigned int jType = iType; jType < hist[0].size(); ++jType) {
            file << histIntra2_[im][iType][jType][bin] << " ";
          }
        }
        for (unsigned int iType = 0; iType < histMoments_[im].size(); ++iType) {
          for (unsigned int jType = iType; jType < histMoments_[im][0].size(); ++jType) {
            for (unsigned int iMo = 0; iMo < histMoments_[im][0][0].size(); ++iMo) {
              file << histMoments_[im][iType][jType][iMo][bin] << " ";
            }
          }
        }
        file << endl;
      }
    }
  }
}

void AnalyzeScatter::defaultConstruction_() {
  className_.assign("AnalyzeScatter");
  nMomentsCut_ = 0;
  verbose_ = 0;
}

void AnalyzeScatter::initSANS(
  const double dgr,
  const int nMacros
  ) {
  countConf_.resize(nMacros, 0.);
  dgr_ = dgr;
  nbins_ = static_cast<int>((space()->minl()/2-10*DTOL)/dgr_) + 1;
  qbins_ = static_cast<int>(nbins_)+1;
  // qbins_ = static_cast<int>(nbins_*10)+1;
  if (qbins_ > 500) qbins_ = 500;
  iqm_.resize(nMacros, vector<double>(qbins_, 0.));
  const double qmin = log(4*PI/space()->minl());
  double sigmin;
  const vector<double> sig = pair_->sig();
  sigmin = *std::min_element(sig.begin(), sig.end());
  const double qmax = log(6*PI/std::min(sigmin, space()->minBondLength())),
    dq = (qmax - qmin)/static_cast<double>(qbins_);
  qwave_.resize(qbins_);
  Pq_.resize(qbins_, vector<double>(nPartTypes_()));
  for (int q = 0; q < qbins_; ++q) {
    const double qtmp = qmin + dq*q;
    qwave_[q] = exp(qtmp);
    for (int itype = 0; itype < nPartTypes_(); ++itype) {
      double r;
      r = 0.5*pair_->sig()[itype];
      const double qr = r*qwave_[q];
      Pq_[q][itype] = 3.*(sin(qr)-qr*cos(qr))/(qr*qr*qr);
    }
  }

  // initialize histograms
  histInter2_.resize(nMacros, vector<vector<vector<long long> > >(
    nPartTypes_(), vector<vector<long long> >(
      nPartTypes_(), vector<long long>(nbins_))));
  histIntra2_.resize(nMacros, vector<vector<vector<long long> > >(
    nPartTypes_(), vector<vector<long long> >(
      nPartTypes_(), vector<long long>(nbins_))));
  histMoments_.resize(nMacros, vector<vector<vector<vector<double> > > >(
    nPartTypes_(), vector<vector<vector<double> > >(
      nPartTypes_(), vector<vector<double> >(
        nMomentsCut_, vector<double>(nbins_)))));
}

void AnalyzeScatter::update(const int iMacro) {
  ++countConf_[iMacro];
  const vector<double> l = space()->l();
  const vector<int> type = space()->type();
  const vector<int> mol = space()->mol();
  const vector<double> x = space()->x();
  const int natom = space()->natom();
  const double minl = space()->minl();
  const int dimen = space()->dimen();

  double lx = l[0], ly = 0, lz = 0, halflx = lx/2., halfly = 0, halflz = 0,
    dx = 0, dy = 0, dz = 0, xi = 0, yi = 0, zi = 0, r2 = 0;
  if (dimen >= 2) ly = l[1], halfly = ly/2.;
  if (dimen >= 3) lz = l[2], halflz = lz/2.;
  for (int ipart = 0; ipart < natom-1; ++ipart) {
    const int iType = type[ipart];
    const int iMol = mol[ipart];
    xi = x[dimen*ipart];
    if (dimen >= 2) yi = x[dimen*ipart+1];
    if (dimen >= 3) zi = x[dimen*ipart+2];
    for (int jpart = ipart + 1; jpart < natom; ++jpart) {
      const int jMol = mol[jpart];
      if (iMol != jMol) {
        const int jType = type[jpart];
        dx = xi - x[dimen*jpart];
        if (dx >  halflx) dx -= lx;
        if (dx < -halflx) dx += lx;
        if (dimen >= 2) {
          dy = yi - x[dimen*jpart+1];
          if (dy >  halfly) dy -= ly;
          if (dy < -halfly) dy += ly;
        }
        if (dimen >= 3) {
          dz = zi - x[dimen*jpart+2];
          if (dz >  halflz) dz -= lz;
          if (dz < -halflz) dz += lz;
        }
        r2 = dx*dx + dy*dy + dz*dz;
        const double r = sqrt(r2);
        if (r <= 0.5*minl) {
          const int bin = feasstRound(r/dgr_ - 0.5);
          histInter2_[iMacro][iType][jType][bin]++;
          histInter2_[iMacro][jType][iType][bin]++;

          // moments
          double peMoment = 1.;
          for (int iMo = 0; iMo < nMomentsCut_; ++iMo) {
            histMoments_[iMacro][iType][jType][iMo][bin] += peMoment;
            histMoments_[iMacro][jType][iType][iMo][bin] += peMoment;
            peMoment *= pair_->peTot();
          }
        }
      }
    }
  }

  // intramolecular term
  const vector<int> mol2part = space()->mol2part();
  for (int iMol = 0; iMol < space()->nMol(); ++iMol) {
    for (int ipart = mol2part[iMol]; ipart < mol2part[iMol+1]-1; ++ipart) {
      const int iType = type[ipart];
      xi = x[dimen*ipart];
      if (dimen >= 2) yi = x[dimen*ipart+1];
      if (dimen >= 3) zi = x[dimen*ipart+2];
      for (int jpart = ipart+1; jpart < mol2part[iMol+1]; ++jpart) {
        const int jType = type[jpart];
        dx = xi - x[dimen*jpart];
        if (dx >  halflx) dx -= lx;
        if (dx < -halflx) dx += lx;
        if (dimen >= 2) {
          dy = yi - x[dimen*jpart+1];
          if (dy >  halfly) dy -= ly;
          if (dy < -halfly) dy += ly;
        }
        if (dimen >= 3) {
          dz = zi - x[dimen*jpart+2];
          if (dz >  halflz) dz -= lz;
          if (dz < -halflz) dz += lz;
        }
        r2 = dx*dx + dy*dy + dz*dz;
        const double r = sqrt(r2);
        if (r <= 0.5*minl) {
          const int bin = feasstRound(r/dgr_ - 0.5);
          histIntra2_[iMacro][iType][jType][bin]++;
          histIntra2_[iMacro][jType][iType][bin]++;
        }
      }
    }
  }
}

void AnalyzeScatter::computeSANS(const int iMacro,
                       const int nMol
  ) {
  iq_.resize(qbins_);
  iqIntra_.resize(qbins_);
  iqInter_.resize(qbins_);
  for (int k = 0; k < qbins_; ++k) {
    const double q = qwave_[k];

    // compute Pq2 for one molecule
    double Pq2 = 0;
    for (int iType = 0; iType < nPartTypes_(); ++iType) {
      Pq2 += pow(Pq_[k][iType], 2);
    }

    // compute the various parts of the intensity by Fourier transform
    double Iq1 = 0, Iq2 = 0, Iq3 = 0;
    for (int bin = 0; bin < nbins_; ++bin) {
      const double r = dgr_*(bin + 0.5),
        qr = q*r,
        sqrinvqr = sin(qr)/qr,
        rmin = dgr_*bin,
        rmax = dgr_*(bin + 1),
        dv = 4./3.*PI*(pow(rmax, 3) - pow(rmin, 3))/space()->vol();
      double nid;
      nid = dv*(nMol - 1);
      for (int iType = 0; iType < nPartTypes_(); ++iType) {
        for (int jType = 0; jType < nPartTypes_(); ++jType) {
          double normFacIn;
          normFacIn = static_cast<double>(nMol);
          Iq1 += Pq_[k][iType]*Pq_[k][jType]*
            histIntra2_[iMacro][iType][jType][bin]*sqrinvqr/normFacIn
              /static_cast<double>(countConf_[iMacro]);
          Iq2 += Pq_[k][iType]*Pq_[k][jType]
            *histInter2_[iMacro][iType][jType][bin]*sqrinvqr/normFacIn
              /static_cast<double>(countConf_[iMacro]);
          Iq3 += Pq_[k][iType]*Pq_[k][jType]*nid*sqrinvqr;
        }
      }
    }
    double normFac = pow(nPartTypes_(), 2);
    Iq1 /= normFac; Iq2 /= normFac; Iq3 /= normFac; Pq2 /= normFac;
    iq_[k] = Pq2 + Iq1 + Iq2 - Iq3;
    iqIntra_[k] = Pq2 + Iq1;
    iqInter_[k] = Iq2 - Iq3;
  }
}

void AnalyzeScatter::write() {
  computeSANS();
  stringstream ss;
  ss << fileName_ << ".txt";
  printer_(ss.str());
}

void AnalyzeScatter::printer_(const string fileName, CriteriaWLTMMC *c,
  const int iMacro) {
  // initialize output
  std::ofstream file(fileName.c_str());
  stringstream ss;
  ss << "#q Iq Intra Inter" << endl;
  if (fileName.empty()) {
    cout << ss.str();
  } else {
    file << ss.str();
  }

  // print
  for (int k = 0; k < qbins_; ++k) {
    ss.str("");
    ss << qwave_[k] << " " << iq_[k] << " " << iqIntra_[k] << " "
       << iqInter_[k] << " " << iq_[k]/iqIntra_[k] << endl;
    if (fileName.empty()) {
      cout << ss.str();
    } else {
      file << ss.str();
    }
  }

  // initialize histogram output
  ss.str("");
  ss << fileName << "h";
  // std::ofstream file2(ss.str().c_str());
  writeRestart(ss.str().c_str(), iMacro);

  // obtain the number of molecules
  int nMol = -1;
  if (c != NULL) {
    if (c->mType().compare("nmol") == 0) {
      nMol = c->bin2m(iMacro);
    } else {
      nMol = space()->nMol();
    }
  } else {
    nMol = space()->nMol();
  }

  const int maxBins = static_cast<int>(space()->minl()/2./dgr_);

  // compute and print radial distribution functions
  ss.str("");
  ss << fileName << "g";
  std::ofstream file5(ss.str().c_str());
  ss.str("");
  ss << "#r ";
  vector<vector<vector<long long> > > &hist = histInter2_[iMacro];
  vector<vector<vector<double> > > gr(hist.size(), vector<vector<double> >(
    hist[0].size(), vector<double>(maxBins, 0.)));
  for (int bin = 0; bin < maxBins; ++bin) {
    ss.str("");
    const double r = dgr_*(bin + 0.5),
      rmin = r - 0.5*dgr_,
      rmax = r + 0.5*dgr_,
      dv = 4./3.*PI*(pow(rmax, space()->dimen())-pow(rmin,
        space()->dimen()))/space()->vol();
    ss << r << " ";
    for (unsigned int iType = 0; iType < hist.size(); ++iType) {
      const int niType = nMol*space()->addMolList()[0]->nType()[iType];
      for (unsigned int jType = 0; jType < hist[0].size(); ++jType) {
        int normFac = 0;
        if (iType == jType) normFac = 1;
        const int njType = nMol*space()->addMolList()[0]->nType()[jType];
        gr[iType][jType][bin] = (hist[iType][jType][bin])
          /static_cast<double>( (niType - normFac)*njType*(countConf_[iMacro]) )
          /dv;
        if (iType <= jType) ss << gr[iType][jType][bin] << " ";
      }
    }
    ss << endl;
    if (fileName.empty()) {
      cout << ss.str();
    } else {
      file5 << ss.str();
    }
  }

  // find scaling factor from last 5% of gr
  int nlast = static_cast<int>(0.05*maxBins);
  vector<vector<double> > sgr(hist.size(), vector<double>(hist[0].size(), 0.));
  for (unsigned int iType = 0; iType < hist.size(); ++iType) {
  for (unsigned int jType = 0; jType < hist[0].size(); ++jType) {
    for (int bin = maxBins - nlast; bin < maxBins; ++bin) {
      sgr[iType][jType] += gr[iType][jType][bin]/static_cast<double>(nlast);
    }
  }}

  // output rescaled gr
  ss.str("");
  ss << fileName << "gs";
  std::ofstream file3(ss.str().c_str());
  for (int bin = 0; bin < maxBins; ++bin) {
    ss.str("");
    const double r = dgr_*(bin + 0.5);
    ss << r << " ";
    for (unsigned int iType = 0; iType < hist.size(); ++iType) {
      for (unsigned int jType = iType; jType < hist[0].size(); ++jType) {
        ss << gr[iType][jType][bin]/sgr[iType][jType] << " ";
      }
    }
    ss << endl;
    if (fileName.empty()) {
      cout << ss.str();
    } else {
      file3 << ss.str();
    }
  }

  // fourier transform again
  ss.str("");
  ss << fileName << "f";
  std::ofstream file4(ss.str().c_str());
  for (int k = 0; k < qbins_; ++k) {
    const double q = qwave_[k];
    ss.str("");
    ss << q << " ";

    // compute Pq2 for one molecule
    double Pq2 = 0;
    for (int iType = 0; iType < nPartTypes_(); ++iType) {
      Pq2 += pow(Pq_[k][iType], 2);
    }

    // compute the various parts of the intensity by Fourier transform
    double Iq1 = 0., Iq2 = 0., Iq2shift = 0., Iq2scale = 0., Iq3 = 0.;
    for (int bin = 0; bin < maxBins; ++bin) {
      const double r = dgr_*(bin + 0.5),
        qr = q*r,
        sqrinvqr = sin(qr)/qr,
        rmin = dgr_*bin,
        rmax = dgr_*(bin + 1),
        dv = 4./3.*PI*(pow(rmax, 3) - pow(rmin, 3))/space()->vol();
      double nid;
      nid = dv*(nMol - 1);
      for (int iType = 0; iType < nPartTypes_(); ++iType) {
        for (int jType = 0; jType < nPartTypes_(); ++jType) {
          double normFacIn;
          normFacIn = static_cast<double>(nMol);
          Iq1 += Pq_[k][iType]*Pq_[k][jType]*
            histIntra2_[iMacro][iType][jType][bin]*sqrinvqr/normFacIn
            /static_cast<double>(countConf_[iMacro]);
          // Iq2 += Pq_[k][iType]*Pq_[k][jType]*
          //  histInter2_[iMacro][iType][jType][bin]*sqrinvqr/normFacIn
          //  /static_cast<double>(countConf_[iMacro]);
          // Iq2 += Pq_[k][iType]*Pq_[k][jType]*
          //  histInter2_[iMacro][iType][jType][bin]*sqrinvqr/normFacIn
          //  /static_cast<double>(countConf_[iMacro])/sgr[iType][jType];
          // Iq3 += Pq_[k][iType]*Pq_[k][jType]*nid*sqrinvqr;
          if (nMol > 1) {
            Iq2 += Pq_[k][iType]*Pq_[k][jType]* (gr[iType][jType][bin]-1.)
                   *sqrinvqr*nid;
            Iq2shift += Pq_[k][iType]*Pq_[k][jType]* (gr[iType][jType][bin]
              -sgr[iType][jType])*sqrinvqr*nid;
            Iq2scale += Pq_[k][iType]*Pq_[k][jType]* (gr[iType][jType][bin]
              /sgr[iType][jType]-1.)*sqrinvqr*nid;
          }
        }
      }
    }
    double normFac = pow(nPartTypes_(), 2);
    // if (nMol == 0) normFac = 1;
    Iq1 /= normFac; Iq2 /= normFac; Iq3 /= normFac; Pq2 /= normFac;
    Iq2shift /= normFac; Iq2scale /= normFac;
    // iq_[k] = Pq2 + Iq1 + Iq2 - Iq3;
    //
    // iqIntra_[k] = Pq2 + Iq1;
    // iqInter_[k] = Iq2 - Iq3;
    ss << Pq2 + Iq1 + Iq2 - Iq3 << " " << Pq2 + Iq1 << " " << Iq2 - Iq3
       << " " << Pq2 + Iq1 + Iq2shift - Iq3 << " "
       << Pq2 + Iq1 + Iq2scale - Iq3 << endl;
    iqm_[iMacro][k] = Pq2 + Iq1 + Iq2scale - Iq3;
    if (fileName.empty()) {
      cout << ss.str();
    } else {
      file4 << ss.str();
    }
  }
}

void AnalyzeScatter::write(CriteriaWLTMMC *c) {
  for (int iMacro = 0; iMacro < c->nBin(); ++iMacro) {
    // obtain the number of molecules
    int nMol = -1;
    if (c->mType().compare("nmol") == 0) {
      nMol = c->bin2m(iMacro);
    } else {
      nMol = space()->nMol();
    }

    // obtain the scattering intensity
    computeSANS(iMacro, nMol);

    stringstream ss;
    ss << fileName_ << "i" << iMacro << ".txt";
    printer_(ss.str(), c, iMacro);
  }
}

int AnalyzeScatter::nPartTypes_() {
  return space()->nParticleTypes();
}

shared_ptr<AnalyzeScatter> makeAnalyzeScatter(Pair *pair, const argtype &args) {
  return make_shared<AnalyzeScatter>(pair, args);
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_



