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
 * The abbreviated user interface allows one-line instantiation of objects
 * such as Monte Carlo trials.
 */

#ifndef UI_ABBREVIATED_H_
#define UI_ABBREVIATED_H_

#include <memory>
#include "./mc.h"
#include "./mc_wltmmc.h"

namespace feasst {

class AbbreviatedUI {};  // here is where we trick ../tools/makeFactory.sh.

/// Initialize TrialAdd and TrialDelete in mc, with equal weights.
void insertDeleteTrial(MC *mc, const char* moltype);

/// Initialize TrialAdd and TrialDelete in mc, with equal weights.
void insertDeleteTrial(shared_ptr<MC> mc, const char* moltype);

}  // namespace feasst

#endif  // UI_ABBREVIATED_H_
