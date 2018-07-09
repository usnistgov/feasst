/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef SRC_BOND_H_
#define SRC_BOND_H_

#include "atom.h"

namespace feasst {

class Space;

/**
 * List first bond per atom/site.
 */
class Bond : public Atom {
 public:
  /// Constructor
  Bond();

 protected:
  // Set value
  void setVal_(const Space &space, const int iatom);
};

}  // namespace feasst

#endif  // SRC_BOND_H_
