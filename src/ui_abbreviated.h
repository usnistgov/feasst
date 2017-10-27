/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

/**
 * \file
 *
 * \brief This is the user interface class
 *
 * The abbreviated user interface allows one-line instantiation of trials for
 * Monte Carlo. Note that while shared_ptr are encourged, one may also use raw
 * pointers to create the Trials; however, this could lead to memory leaks
 * and other issues. Shared pointers are not as simple to use for the python
 * interface.
 *
 */

#ifndef UI_ABBREVIATED_H_
#define UI_ABBREVIATED_H_

#include <memory>
#include "./mc.h"
#include "./mc_wltmmc.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

class AbbreviatedUI {};  // here is where we trick ../tools/makeFactory.sh.

void insertDeleteTrial(MC *mc, const char* moltype);
void insertDeleteTrial(shared_ptr<MC> mc, const char* moltype);

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // UI_ABBREVIATED_H_
