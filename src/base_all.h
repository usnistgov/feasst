#ifndef BASE_ALL_H_
#define BASE_ALL_H_

#include "./base_random.h"

namespace feasst {

/**
 * This class is a container for all other classes to inherit.
 */
class BaseAll : public BaseRandom {
 public:
  BaseAll();
  virtual ~BaseAll() {}

 protected:
};

}  // namespace feasst

#endif  // BASE_ALL_H_

