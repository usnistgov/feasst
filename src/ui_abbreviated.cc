/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef UI_ABBREVIATED_H_
#define UI_ABBREVIATED_H_

#include "./ui_abbreviated.h"
#include "./trial_add.h"
#include "./trial_delete.h"
#include "./mc.h"

//using namespace feasst;
namespace feasst {

void insertDeleteTrial(MC *mc, const char* moltype) {
  addTrial(mc, moltype);
  deleteTrial(mc, moltype);
}
void insertDeleteTrial(shared_ptr<MC> mc, const char* moltype) {
  insertDeleteTrial(mc.get(), moltype);
}

}  // namespace feasst

#endif  // UI_ABBREVIATED_H_
