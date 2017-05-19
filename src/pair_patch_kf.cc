#include "./pair_patch_kf.h"
#include "./mins.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

PairPatchKF::PairPatchKF(Space *space,
  const double rCut,        //!< interaction cut-off distance (square well)
  const double patchAngle    //!< solid half-angle of patch (see reference)
  )
  : Pair(space, rCut),
    patchAngle_(patchAngle),
    cpa_(cos(patchAngle/180*PI)) {
  className_.assign("PairPatchKF");
  mirrorPatch(0);
  initAtomCut(0);
}

PairPatchKF::~PairPatchKF() {
}

/**
 * pair-wise force calculation of entire configuration
 */
int PairPatchKF::initEnergy() {
  // shorthand for read-only space variables
  const int nMol = space_->nMol();
  const vector<double> &l = space_->l();
  const vector<double> &x = space_->x();
  const vector<int> &mol2part = space_->mol2part();

  // zero accumulators: potential energy, force, and virial
  std::fill(pe_.begin(), pe_.end(), 0.);
  fill(0., f_);
  fill(0., vr_);
  peSR_ = 0;

  double pePart = 0;
  const double sigSq = pow(sig_[0], 2);

  // loop through pairs of molecules
  for (int iMol = 0; iMol < nMol - 1; ++iMol) {
    const int ipart = mol2part[iMol];
    for (int jMol = iMol + 1; jMol < nMol; ++jMol) {
      const int jpart = mol2part[jMol];
      // separation vector, xij with periodic boundary conditions
      double r, r2 = 0;
      vector<double> xij(dimen_);
      for (int i = 0; i < dimen_; ++i) {
        r = x[dimen_*ipart+i] - x[dimen_*jpart+i];
        if (r >  0.5 * l[i]) r -= l[i];
        if (r < -0.5 * l[i]) r += l[i];
        xij[i] = r;
        r2 += r*r;
      }

      // no interaction beyond cut-off distance
      if (r2 < rCutSq_) {
        // hard sphere
        if (r2 < sigSq) {
          pePart = std::numeric_limits<double>::max()/1e10;

        // orientational square well
        } else {
          // loop through pairs of patches to see how many patch-patch
          // configurations are present
          const double r = sqrt(r2);
          pePart = 0;
          for (int isite = ipart + 1; isite < mol2part[iMol+1]; ++isite) {
            for (int jsite = jpart + 1; jsite < mol2part[jMol+1]; ++jsite) {
              // cosangle between xij and ipatch
              double cosa = 0;
              for (int dim = 0; dim < dimen_; ++dim) {
                cosa -= xij[dim]*(x[dimen_*isite+dim] - x[dimen_*ipart+dim]);
              }
              cosa /= r;
              if ( (cosa >= cpa_) || ( (mirrorPatch_) && (cosa <= -cpa_) ) ) {
                // cosangle between xij and jpatch
                cosa = 0;
                for (int dim = 0; dim < dimen_; ++dim) {
                  cosa += xij[dim]*(x[dimen_*jsite+dim] - x[dimen_*jpart+dim]);
                }
                cosa /= r;
                if ( (cosa >= cpa_) || ( (mirrorPatch_) && (cosa <= -cpa_) ) ) {
                  pePart -= 1.;
                }
              }
            }
          }
        }
        peSR_ += pePart;
      }
    }
  }
  return 0;
}

/**
 * Lennard-Jones potential energy contribution due to particles
 */
double PairPatchKF::multiPartEner(
  const vector<int> mpart,    //!< particles to calculate energy interactions
  const int flag     //!< place holder for other pair styles
  ) {
  if (flag == 0) {}  // remove unused parameter warning

  // zero potential energy contribution of particle ipart
  peSRone_ = 0;

  /// obtain neighList with cellList
  if (space_->cellType() == 1) {
    space_->buildNeighListCell(space_->mol()[mpart.front()]);
  }

  return multiPartEnerNeigh(mpart);
}

/**
 * Lennard-Jones potential energy contribution due to particles
 */
double PairPatchKF::multiPartEnerNeigh(
  const vector<int> mpart    //!< particles to calculate energy interactions
  ) {
  // shorthand for read-only space variables
  const vector<double> &x = space_->x();
  const vector<double> &l = space_->l();
  const vector<int> &mol2part = space_->mol2part();

  // declare variables for optimization
  int jpart, isite, jsite;
  double cosa, r2, dx, dy, dz;
  const double lx = l[0], ly = l[1], lz = l[2],
    halflx = lx/2., halfly = ly/2., halflz = lz/2.;
  const int ipart = mpart[0], iMol = space_->mol()[ipart];
  const double xi = x[dimen_*ipart+0];
  const double yi = x[dimen_*ipart+1];
  const double zi = x[dimen_*ipart+2];
  const double sigSq = pow(sig_[0], 2), cpaSq = cpa_*cpa_;

  // avoid the sqrt operation for patches of solid angle <= 90
  //  first, compute angle*|rij| and check its sign.
  //  if 90 degree patch, its a patch of >= 0
  //  if <90 degree patch, compute angle = (angle*|rij|)^2/rij2
  ASSERT(patchAngle_ <= 90,
    "Optimization built for patches of solid angle <= 90");

  if (neighOn_) {
    neighOne_.clear();
    neighOne_.resize(static_cast<int>(mpart.size()), vector<int>());
  }

  // loop through molecules
  const vector<int> &neigh = space_->neighListChosen();
  for (unsigned int ineigh = 0; ineigh < neigh.size(); ++ineigh) {
    const int jMol = neigh[ineigh];
    if (iMol != jMol) {
      jpart = mol2part[jMol];
      // separation distance with periodic boundary conditions
      dx = xi - x[dimen_*jpart];
      dy = yi - x[dimen_*jpart+1];
      dz = zi - x[dimen_*jpart+2];
      if (dx >  halflx) dx -= lx;
      if (dx < -halflx) dx += lx;
      if (dy >  halfly) dy -= ly;
      if (dy < -halfly) dy += ly;
      if (dz >  halflz) dz -= lz;
      if (dz < -halflz) dz += lz;
      r2 = dx*dx + dy*dy + dz*dz;

      // no interaction beyond cut-off distance
      if (r2 < rCutSq_) {
        // store new neighbor list
        if (neighOn_) {
          if ( (r2 < neighAboveSq_) && (r2 > neighBelowSq_) ) {
            neighOne_[0].push_back(jMol);
          }
        }

        // hard sphere
        if (r2 < sigSq) {
          peSRone_ += std::numeric_limits<double>::max()/1e10;

        // orientational square well
        } else {
          // loop through pairs of patches
          for (isite = ipart + 1; isite < mol2part[iMol+1]; ++isite) {
            for (jsite = jpart + 1; jsite < mol2part[jMol+1]; ++jsite) {
//              const double r = sqrt(r2);
//              cosa = ( dx*(xi - x[dimen_*isite]  ) +
//                       dy*(yi - x[dimen_*isite+1]) +
//                       dz*(zi - x[dimen_*isite+2]) ) / r;
//              if ( (cosa >=  cpa_) || ( (mirrorPatch_) && (cosa <= -cpa_) ) ){
//                cosa = ( dx*(x[dimen_*jsite]   - x[dimen_*jpart]) +
//                         dy*(x[dimen_*jsite+1] - x[dimen_*jpart+1]) +
//                         dz*(x[dimen_*jsite+2] - x[dimen_*jpart+2]) ) / r;
//              }
              cosa = dx*(xi - x[dimen_*isite])   +
                     dy*(yi - x[dimen_*isite+1]) +
                     dz*(zi - x[dimen_*isite+2]);
              if ( (mirrorPatch_) || (cosa >= 0) ) {
                cosa = cosa*cosa/r2;
                if (cosa >=  cpaSq) {
                  cosa = dx*(x[dimen_*jsite]   - x[dimen_*jpart]) +
                         dy*(x[dimen_*jsite+1] - x[dimen_*jpart+1]) +
                         dz*(x[dimen_*jsite+2] - x[dimen_*jpart+2]);
                  if ( (mirrorPatch_) || (cosa >= 0) ) {
                    cosa = cosa*cosa/r2;
                    if (cosa >=  cpaSq) peSRone_ -= 1.;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
//  cout << peSRone_ << endl;
  return peSRone_;
}

/**
 * stores, restores or updates variables to avoid recompute of entire
 * configuration after every change
 */
void PairPatchKF::update(
  const vector<int> mpart,    //!< particles involved in move
  const int flag,         //!< type of move
  const char* uptype    //!< description of update type
  ) {
  if (neighOn_) {
    updateBase(mpart, flag, uptype, neigh_, neighOne_, neighOneOld_);
//  updateBase(mpart, flag, uptype, neighCut_, neighCutOne_, neighCutOneOld_);
  }
  std::string uptypestr(uptype);

  if (uptypestr.compare("store") == 0) {
    if (flag == 0 || flag == 2 || flag == 3) {
      deSR_ = peSRone_;
    }
  }

  if (uptypestr.compare("update") == 0) {
    if (flag == 0) {
      peSR_ += peSRone_ - deSR_;
    }
    if (flag == 2) {
      peSR_ -= deSR_;
    }
    if (flag == 3) {
      peSR_ += deSR_;
    }
  }
}

/**
 * write xyz for visualization
 *  use when mirrorPatch
 */
int PairPatchKF::printxyz(const char* fileName,   //!< file with configuration
  const int initFlag,    //!< open if flag is 1, append if flag is 0
  const std::string comment) {
  ASSERT(dimen_ == 3, "printxyz assumes three dimensions");

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
  const double r1 = sig_[0]/2., eps = r1/100.,
  r2 = (2*r1*(1-cpa_)*(r1+eps)+eps*eps)/(2*eps+2*r1*(1-cpa_)),
  a = r1 + eps - r2;

  int natom = space_->natom();
  if (mirrorPatch_) {
    natom += space_->natom() - space_->nMol();
  }
  fprintf(xyzFile, "%d\n1 %s\n", natom, comment.c_str());
  if (xyzFile != NULL) {
    for (int iMol = 0; iMol < space_->nMol(); ++iMol) {
      const int ipart = space_->mol2part()[iMol];
      for (int isite = ipart; isite < space_->mol2part()[iMol+1]; ++isite) {
        if (isite != ipart) {
          fprintf(xyzFile, "N ");
          for (int i = 0; i < dimen_; ++i) {
            fprintf(xyzFile, "%f ", space_->x(ipart, i)
              + a*(space_->x(isite, i) - space_->x(ipart, i)));
          }
          fprintf(xyzFile, "\n");
          if (mirrorPatch_) {
            fprintf(xyzFile, "N ");
            for (int i = 0; i < dimen_; ++i) {
              fprintf(xyzFile, "%f ", space_->x(ipart, i)
                + a*(-space_->x(isite, i) + space_->x(ipart, i)));
            }
            fprintf(xyzFile, "\n");
          }
        } else {
          fprintf(xyzFile, "O ");
          for (int i = 0; i < dimen_; ++i) {
            fprintf(xyzFile, "%f ", space_->x(isite, i));
          }
          fprintf(xyzFile, "\n");
        }
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

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_






