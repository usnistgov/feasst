/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./base.h"
#include <stdio.h>
#include <iostream>
#include <numeric>

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

Base::Base() : verbose_(0) {
  className_.assign("Base");
  stringstream ss;
  ss << TOSTRING(FEASST_SRC_) << "/..";
  install_dir_.assign(ss.str().c_str());
}

void Base::reconstruct() {
  reconstructDerived_();
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

