/**
 * \file
 *
 * \brief
 */

#include "./base.h"
#include <stdio.h>
#include <iostream>
#include <numeric>

/**
 * Constructor
 */
Base::Base() : verbose_(0) {
  className_.assign("Base");
  stringstream ss;
  ss << TOSTRING(FEASST_SRC_) << "/..";
  install_dir_.assign(ss.str().c_str());
}

/**
 * reset object pointers
 */
void Base::reconstruct() {
  reconstructDerived();
}

