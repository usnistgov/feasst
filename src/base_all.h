/**
 * \file
 *
 * \brief interface for all classes to inherit
 *
 */

#ifndef BASEALL_H_
#define BASEALL_H_

#include "./base_math.h"

namespace feasst {

class BaseAll : public BaseMath {
 public:
  BaseAll();
  virtual ~BaseAll() {}

 protected:
};

}  // namespace feasst

#endif  // BASEALL_H_

