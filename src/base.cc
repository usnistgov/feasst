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

