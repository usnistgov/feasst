/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./bond.h"
#include "./space.h"

namespace feasst {

Bond::Bond() : Atom() {}

void Bond::setVal_(const Space &space, const int iatom) {
  vector<vector<int> > bonds = space.listBonds(iatom);
  if (bonds.size() <= 0) {
    return;
  }
  const int iMol = space.mol()[iatom];
  const int firstAtom = space.mol2part()[iMol];
  int bonded = bonds[0][1];
  if (bonded == iatom - firstAtom) bonded = bonds[0][2];
  intVal_[iatom] = bonded;
}

}  // namespace feasst
