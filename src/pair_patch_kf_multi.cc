/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./pair_patch_kf_multi.h"
#include <math.h>

// ADDED REINHARDT FUNCTION
namespace feasst {

PairPatchKFMulti::PairPatchKFMulti(Space *space, const argtype &args)
  : PairPatchKF(space, args) {
  initAtomCut(1);
  cpaSqi_.resize(2);
  cpaSqi_[0] = cpa_;
  cpaSqi_[1] = cpa_;
  bond_ = make_shared<Bond>();
  space->initAtom(bond_);
}

void PairPatchKFMulti::initPatchAngle(const double angle) {
  PairPatchKF::initPatchAngle(angle);
  initIJ();
}

// templated from pairLoopSite_ in pair.cc
double PairPatchKFMulti::multiPartEnerNeigh(
  const vector<int> siteList) {
  // shorthand for read-only space variables
  const vector<int> &type = space_->type();
  const vector<int> &mol = space_->mol();
  const vector<double> &x = space_->x();
  const vector<double> &boxLength = space_->boxLength();

  // declare variables for optimization
  double dx, dy, dz = 0., zi = 0.;
  double cosi, cosj;

  // PBC optimization variables
  const double lx = boxLength[0];
  const double ly = boxLength[1];
  double lz = 0.;
  if (dimen_ >= 3) {
    lz = boxLength[2];
  }
  const double xyTilt = space_->xyTilt();
  const double xzTilt = space_->xzTilt();
  const double yzTilt = space_->yzTilt();
  const double halflx = lx/2., halfly = ly/2., halflz = lz/2.;
  const int xpbc = space_->periodic(0),
            ypbc = space_->periodic(1);
  int zpbc = 0;
  if (dimen_ >= 3) {
    zpbc = space_->periodic(2);
  }

  // avoid the sqrt operation for patches of solid angle <= 90
  //  first, compute angle*|rij| and check its sign.
  //  if 90 degree patch, its a patch of >= 0
  //  if <90 degree patch, compute angle = (angle*|rij|)^2/rij2
  ASSERT(patchAngle_ <= 90,
    "Optimization built for patches of solid angle <= 90");

  // initialize neigh, neighCut and peMap
  initNeighCutPEMap(siteList);

  // loop through all sites
  for (unsigned int ii = 0; ii < siteList.size(); ++ii) {
    const int iSite = siteList[ii];
    const int itype = type[iSite];
    if (fabs(eps_[itype]) > DTOL) {
      const int iMol = space_->mol()[iSite];
      const int iFirstSite = space_->mol2part()[iMol];
      const double isig = sig_[itype];
      const double xi = x[dimen_*iSite],
                   yi = x[dimen_*iSite+1];
      if (dimen_ >= 3) {
         zi = x[dimen_*iSite+2];
      }

      // obtain neighList with cellList
      if (useCellForSite_(itype) ) {
        space_->buildNeighListCellAtomCut(iSite);
      } else {
        space_->initAtomCut(1);   // set neighListChosen to all atoms
      }
      const vector<int> &neigh = space_->neighListChosen();

      // loop neighboring sites
      for (unsigned int ineigh = 0; ineigh < neigh.size(); ++ineigh) {
        const int jSite = neigh[ineigh];
        const int jMol = mol[jSite];
        const int jFirstSite = space_->mol2part()[jMol];
        const int jtype = type[jSite];
        if ( (iMol != jMol) && (fabs(eps_[jtype]) > DTOL) ) {
          const double jsig = sig_[jtype];
          // compute patch-patch interaction if both sigmas are zero
          const bool patch = (fabs(isig) < DTOL) && (fabs(jsig) < DTOL);
          // compute hard-sphere interaction if both sigmas are nonzero
          const bool hard = (fabs(isig) > DTOL) && (fabs(jsig) > DTOL);
          if (hard) {
            // separation distance with periodic boundary conditions
            dx = xi - x[dimen_*jSite];
            dy = yi - x[dimen_*jSite + 1];
            if (dimen_ >= 3) {
              dz = zi - x[dimen_*jSite + 2];
            }
            // optimized macro for PBC
            TRICLINIC_PBC(dx, dy, dz, lx, ly, lz, halflx, halfly, halflz,
                          xyTilt, xzTilt, yzTilt, xpbc, ypbc, zpbc);
            const double r2 = dx*dx + dy*dy + dz*dz;
            const double rCut = rCutij_[itype][jtype];
            if (r2 < rCut*rCut) {
              peSRone_ += NUM_INF;
              storeNeighCutPEMap(jSite, 0);
            }
          } else if (patch) {
            const int iSiteBase = bond_->intVal(iSite) + iFirstSite,
                      jSiteBase = bond_->intVal(jSite) + jFirstSite;
            const double xib = x[dimen_*iSiteBase],
                         yib = x[dimen_*iSiteBase + 1];
            // separation distance with periodic boundary conditions
            dx = x[dimen_*iSiteBase    ] - x[dimen_*jSiteBase];
            dy = x[dimen_*iSiteBase + 1] - x[dimen_*jSiteBase + 1];
            if (dimen_ >= 3) {
              dz = x[dimen_*iSiteBase + 2] - x[dimen_*jSiteBase + 2];
            }
            // optimized macro for PBC
            TRICLINIC_PBC(dx, dy, dz, lx, ly, lz, halflx, halfly, halflz,
                          xyTilt, xzTilt, yzTilt, xpbc, ypbc, zpbc);
            const double r2 = dx*dx + dy*dy + dz*dz;
            const double rCut = rCutij_[itype][jtype];
            if (r2 < rCut*rCut) {
              cosi = dx*(xib - x[dimen_*iSite])   +
                     dy*(yib - x[dimen_*iSite + 1]);
              if (dimen_ >= 3) {
                cosi += dz*(x[dimen_*iSiteBase + 2] - x[dimen_*iSite + 2]);
              }
              if ( (mirrorPatch_) || (cosi >= 0) ) {
                cosi = cosi*cosi/r2;
                if (cosi >=  cpaSqi_[itype]) {
                  const double xjb = x[dimen_*jSiteBase],
                               yjb = x[dimen_*jSiteBase + 1];
                  cosj = dx*(x[dimen_*jSite] - xjb) +
                         dy*(x[dimen_*jSite + 1] - yjb);
                  if (dimen_ >= 3) {
                    cosj += dz*(x[dimen_*jSite + 2] - x[dimen_*jSiteBase + 2]);
                  }
                  if ( (mirrorPatch_) || (cosj >= 0) ) {
                    cosj = cosj*cosj/r2;
                    if (cosj >=  cpaSqi_[jtype]) {
                      // ANN: add function of cosi,cosj
                        if (patchType == 0) {
                          peSRone_ -= epsij_[itype][jtype];
                        } else if (patchType == 1) {
                          // Reinhardt patch with alpha=0
                          peSRone_ -= (epsij_[itype][jtype])*(1./(1.+exp(-10.*cosi*cosi)))*(1./(1.+exp(-10.*cosj*cosj)));
                        } else {
                          ASSERT(0, "Unrecognized patch type");
                        }
                      storeNeighCutPEMap(jSite, 0);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // if using PEMap, peSRone_ was used to map individual potential energies
  // and the total is put in peSRoneAlt_
  if (peMapOn_ == 1) peSRone_ = peSRoneAlt_;

//  cout << peSRone_ << endl;
  return peSRone_;
}

int PairPatchKFMulti::printXYZ(const char* fileName,
  const int initFlag,
  const std::string comment) {
  //ASSERT(dimen_ == 3, "printxyz assumes three dimensions");

  FILE * xyzFile = NULL;
  if (initFlag == 1) {
    fileBackUp(fileName);
    xyzFile = fopen(fileName, "w");
  } else if (initFlag == 0) {
    xyzFile = fopen(fileName, "a");
  } else {
    ASSERT(0, "unrecognized initFlag");
  }

  // make a patch by inscribing sphere of radius r2 inside bead, radius r1=sig
  //  distance between center of sphere and bead is a
  //  for given patch angle, use law of cosines to derive r2
  //  constraint r1 + eps = r2 + a, where eps is small, to avoid clipping
  const double r1 = sig_[0]/2., eps = r1/20.,
  r2 = (2*r1*(1-cpa_)*(r1+eps)+eps*eps)/(2*eps+2*r1*(1-cpa_)),
  a = r1 + eps - r2;

  int natom = 0;
  for (int iSite = 0; iSite < space_->natom(); ++iSite) {
    const int iType = space()->type()[iSite];
    if (fabs(eps_[iType]) > DTOL) {
      ++natom;
    }
  }

  const double lx = space()->boxLength(0);
  double ly = 0, lz = 0, xyt = 0, xzt = 0, yzt = 0;
  if (space()->dimen() > 1) {
    ly = space()->boxLength(1);
    xyt = space()->xyTilt();
  }
  if (space()->dimen() > 2) {
    lz = space()->boxLength(2);
    xzt = space()->xzTilt();
    yzt = space()->yzTilt();
  }
  fprintf(xyzFile, "%d\n1 %f %f %f %f %f %f %s\n", natom, lx, ly, lz,
          xyt, xzt, yzt, comment.c_str());
  if (xyzFile != NULL) {
    for (int iSite = 0; iSite < space_->natom(); ++iSite) {
      const int iType = space()->type()[iSite];
      if (fabs(eps_[iType]) > DTOL) {
        if (fabs(sig_[iType]) < DTOL) {
          const int iMol = space_->mol()[iSite];
          const int iFirstSite = space_->mol2part()[iMol];
          const int iSiteBase = bond_->intVal(iSite) + iFirstSite;
          // draw patches
          fprintf(xyzFile, "N ");
          for (int i = 0; i < dimen_; ++i) {
            fprintf(xyzFile, "%f ", space_->x(iSiteBase, i)
              + a*(space_->x(iSite, i) - space_->x(iSiteBase, i)));
          }
          ASSERT(!mirrorPatch_, "Not implemented for mirror patch.");
        } else {
          fprintf(xyzFile, "O ");
          for (int i = 0; i < dimen_; ++i) {
            fprintf(xyzFile, "%f ", space_->x(iSite, i));
          }
        }
        if (dimen_ == 2) {
          fprintf(xyzFile, " 0");
        }
        fprintf(xyzFile, "\n");
      }
    }
  }
  fclose(xyzFile);

  // write vmd script to visualize in fileName appended with .vmd
  std::stringstream vmdfnamess;
  vmdfnamess << fileName << ".vmd";
  std::ofstream vmdf(vmdfnamess.str().c_str());
  vmdf << "display projection Orthographic" << endl
    << "color Display Background white" << endl
    << "axes location Off" << endl
    << "topo readvarxyz " << trim("/", fileName) << endl
    << "mol modstyle 0 0 VDW 1.0000000 120.000000" << endl
    << "set sel [atomselect top \"name N\"]" << endl
    << "$sel set radius " << r2 << endl
    << "$sel set mass 1" << endl
    << "set sel [atomselect top \"name O\"]" << endl
    << "$sel set radius " << r1 << endl
    << "$sel set mass 1" << endl;

  return 0;
}

double PairPatchKFMulti::allPartEnerForce(const int flag) {
  peSRone_ = 0.;
  if (flag == 0) {
    peSRone_ = peTot();
    return peSRone_;
  }

  // initalize contact map
  contact_.clear();
  contactpbc_.clear();
  contact_.resize(space_->nMol());
  contactpbc_.resize(space_->nMol());
  vector<double> cpbctmplt(space_->dimen());

  // shorthand for read-only space variables
  const vector<int> type = space_->type();
  const vector<int> &mol = space_->mol();
  const vector<double> &x = space_->x();
  const vector<double> &boxLength = space_->boxLength();

  // declare variables for optimization
  double dx, dy, dz = 0., zi = 0.;
  double cosi, cosj;

  // PBC optimization variables
  const double lx = boxLength[0];
  const double ly = boxLength[1];
  double lz = 0.;
  if (dimen_ >= 3) {
    lz = boxLength[2];
  }
  const double xyTilt = space_->xyTilt();
  const double xzTilt = space_->xzTilt();
  const double yzTilt = space_->yzTilt();
  const double halflx = lx/2., halfly = ly/2., halflz = lz/2.;
  const int xpbc = space_->periodic(0),
            ypbc = space_->periodic(1);
  int zpbc = 0;
  if (dimen_ >= 3) {
    zpbc = space_->periodic(2);
  }


  // avoid the sqrt operation for patches of solid angle <= 90
  //  first, compute angle*|rij| and check its sign.
  //  if 90 degree patch, its a patch of >= 0
  //  if <90 degree patch, compute angle = (angle*|rij|)^2/rij2
  ASSERT(patchAngle_ <= 90,
    "Optimization built for patches of solid angle <= 90");

  // loop through all sites
  for (int iSite = 0; iSite < space_->natom() - 1; ++iSite) {
    const int itype = type[iSite];
    if (fabs(eps_[itype]) > DTOL) {
      const int iMol = space_->mol()[iSite];
      const int iFirstSite = space_->mol2part()[iMol];
      const double isig = sig_[itype];
      const double xi = x[dimen_*iSite],
                   yi = x[dimen_*iSite+1];
      if (dimen_ >= 3) {
         zi = x[dimen_*iSite+2];
      }
      for (int jSite = iSite + 1; jSite < space_->natom(); ++jSite) {
        const int jMol = mol[jSite];
        const int jFirstSite = space_->mol2part()[jMol];
        const int jtype = type[jSite];
        if ( (iMol != jMol) && (fabs(eps_[jtype]) > DTOL) ) {
          const double jsig = sig_[jtype];
          // compute patch-patch interaction if both sigmas are zero
          const bool patch = (fabs(isig) < DTOL) && (fabs(jsig) < DTOL);
          // compute hard-sphere interaction if both sigmas are nonzero
          const bool hard = (fabs(isig) > DTOL) && (fabs(jsig) > DTOL);
          if (hard) {
            // separation distance with periodic boundary conditions
            dx = xi - x[dimen_*jSite];
            dy = yi - x[dimen_*jSite + 1];
            if (dimen_ >= 3) {
              dz = zi - x[dimen_*jSite + 2];
            }
            // optimized macro for PBC
            TRICLINIC_PBC(dx, dy, dz, lx, ly, lz, halflx, halfly, halflz,
                          xyTilt, xzTilt, yzTilt, xpbc, ypbc, zpbc);
            const double r2 = dx*dx + dy*dy + dz*dz;
            const double rCut = rCutij_[itype][jtype];
            if (r2 < rCut*rCut) {
              peSRone_ += NUM_INF;
              storeNeighCutPEMap(jSite, 0);
            }
          } else if (patch) {
            const int iSiteBase = bond_->intVal(iSite) + iFirstSite,
                      jSiteBase = bond_->intVal(jSite) + jFirstSite;
            const double xib = x[dimen_*iSiteBase],
                         yib = x[dimen_*iSiteBase + 1];
            // separation distance with periodic boundary conditions
            dx = x[dimen_*iSiteBase    ] - x[dimen_*jSiteBase];
            dy = x[dimen_*iSiteBase + 1] - x[dimen_*jSiteBase + 1];
            if (dimen_ >= 3) {
              dz = x[dimen_*iSiteBase + 2] - x[dimen_*jSiteBase + 2];
            }
            // optimized macro for PBC
            const double dx0 = dx, dy0 = dy, dz0 = dz;
            TRICLINIC_PBC(dx, dy, dz, lx, ly, lz, halflx, halfly, halflz,
                          xyTilt, xzTilt, yzTilt, xpbc, ypbc, zpbc);
            const double r2 = dx*dx + dy*dy + dz*dz;
            const double rCut = rCutij_[itype][jtype];
            if (r2 < rCut*rCut) {
              cosi = dx*(xib - x[dimen_*iSite])   +
                     dy*(yib - x[dimen_*iSite + 1]);
              if (dimen_ >= 3) {
                cosi += dz*(x[dimen_*iSiteBase + 2] - x[dimen_*iSite + 2]);
              }
              if ( (mirrorPatch_) || (cosi >= 0) ) {
                cosi = cosi*cosi/r2;
                if (cosi >=  cpaSqi_[itype]) {
                  const double xjb = x[dimen_*jSiteBase],
                               yjb = x[dimen_*jSiteBase + 1];
                  cosj = dx*(x[dimen_*jSite] - xjb) +
                         dy*(x[dimen_*jSite + 1] - yjb);
                  if (dimen_ >= 3) {
                    cosj += dz*(x[dimen_*jSite + 2] - x[dimen_*jSiteBase + 2]);
                  }
                  if ( (mirrorPatch_) || (cosj >= 0) ) {
                    cosj = cosj*cosj/r2;
                    if (cosj >=  cpaSqi_[jtype]) {
                      const double epsij = epsij_[itype][jtype];
                      // ANN: add function of cosi,cosj
                      if (patchType == 0) {
                        peSRone_ -= epsij;
                      } else if (patchType == 1) {
                        peSRone_ -= epsij*(1./(1.+exp(-10.*cosi*cosi)))*(1./(1.+exp(-10.*cosj*cosj)));
                      } else {
                        ASSERT(0, "Unrecognized patch type");
                      }
                      if (fabs(epsij) > clusterTolerance + DTOL) {
                        storeNeighCutPEMap(jSite, 0);
                        contact_[iMol].push_back(jMol);
                        contact_[jMol].push_back(iMol);
                        cpbctmplt[0] = dx - dx0;
                        cpbctmplt[1] = dy - dy0;
                        if (dimen_ >= 3) {
                          cpbctmplt[2] = dz - dz0;
                        }
                        contactpbc_[iMol].push_back(cpbctmplt);
                        cpbctmplt[0] *= -1;
                        cpbctmplt[1] *= -1;
                        if (dimen_ >= 3) {
                          cpbctmplt[2] *= -1;
                        }
                        contactpbc_[jMol].push_back(cpbctmplt);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return peSRone_;
}

void PairPatchKFMulti::initIJ() {
  cpaSqi_.clear();
  cpaSqi_.resize(space()->nParticleTypes());
  for (int itype = 0; itype < space()->nParticleTypes(); ++itype) {
    cpaSqi_[itype] = cpa_*cpa_;
    const double isig = sig_[itype];
    for (int jtype = itype; jtype < space()->nParticleTypes(); ++jtype) {
      const double jsig = sig_[jtype];
      const double sigij = sigij_[itype][jtype];
      // compute patch-patch interaction if both sigmas are zero
      const bool patch = (fabs(isig) < DTOL) && (fabs(jsig) < DTOL);
      // compute hard-sphere interaction if both sigmas are nonzero
      const bool hard = (fabs(isig) > DTOL) && (fabs(jsig) > DTOL);
      double rc = 0.;
      if (patch) {
        rc = rCut_;
      } else if (hard) {
        rc = sigij;
      }
      // cout << "i " << itype << " " << jtype << " " << rc << endl;
      rCutijset(itype, jtype, rc);
    }
  }
}

void PairPatchKFMulti::initPatchAngleInDegrees(const double angle,
  const int itype) {
  NOTE("A custom angle will be overridden "
    << "during an expanded ensemble patch angle simulation.");
  ASSERT(itype < static_cast<int>(cpaSqi_.size()), "itype(" << itype <<
    ") is too large. "
    << "Did you use initIJ()? Because the patch angle vector is of size "
    << cpaSqi_.size());
  cpaSqi_[itype] = cos(angle/180.*PI);
}

void PairPatchKFMulti::initData(const string fileName) {
  Pair::initData(fileName);
  ASSERT(space()->addMolList().back()->bondList().size() != 0,
    "Bonds are required between patch orienters and patch centers.");

  // check that all patches have a bond and it has unit length
  auto spaceClone = space()->cloneShrPtr();
  ASSERT(spaceClone->nMol() == 0, "initData assumed to be called with no particles");
  PairPatchKFMulti * pair = (*this).clone(spaceClone.get());
  pair->addMol(fileName);
  for (int iSite = 0; iSite < spaceClone->natom(); ++iSite) {
    const int iType = spaceClone->type()[iSite];
    if (fabs(eps_[iType]) > DTOL && fabs(sig_[iType]) < DTOL) {
      const int iSiteBase = pair->bond_->intVal(iSite);
      ASSERT(iSiteBase != pair->bond_->defaultInt(),
        "site(" << iSite <<") is a patch but no bond assigned in order to "
        << "find its center.");
      double r2 = 0.;
      for (int dim = 0; dim < dimen_; ++dim) {
        r2 += pow(spaceClone->x(iSite, dim) - spaceClone->x(iSiteBase, dim), 2);
      }
      const double r = sqrt(r2);
      ASSERT(fabs(r - 1.) < 1e-8,
        "Site(" << iSite << ") is considered a patch unit vector with center "
        << "site(" << iSiteBase << "). The distance between them is "
        << MAX_PRECISION << r << " when it should be unity");
    }
  }
  delete pair;
}

void PairPatchKFMulti::reconstruct(Space* space) {
  bond_ = make_shared<Bond>(*bond_);
  Pair::reconstruct(space);
  space->delPerAtom();
  space->initAtom(bond_);
}

shared_ptr<PairPatchKFMulti> makePairPatchKFMulti(Space *space,
  const argtype &args) {
  return make_shared<PairPatchKFMulti>(space, args);
}

void PairPatchKFMulti::setOrder(const double order) {
  PairPatchKF::setOrder(order);
  if (orderName_.compare("rCut") == 0) {
    rCut_ = order;
    initIJ();
  }
}

}  // namespace feasst
