#ifndef UI_ABBREVIATED_H_
#define UI_ABBREVIATED_H_

#include "./ui_abbreviated.h"
#include "./trial_add.h"
#include "./trial_delete.h"
#include "./mc.h"

//using namespace feasst;
#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

void insertDeleteTrial(MC *mc, const char* moltype) {
  addTrial(mc, moltype);
  deleteTrial(mc, moltype);
}
void insertDeleteTrial(shared_ptr<MC> mc, const char* moltype) {
  insertDeleteTrial(mc.get(), moltype);
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // UI_ABBREVIATED_H_
