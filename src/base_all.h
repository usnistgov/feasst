#ifndef BASE_ALL_H_
#define BASE_ALL_H_

#include "./base_random.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * This class is a container for all other classes to inherit.
 */
class BaseAll : public BaseRandom {
 public:
  BaseAll();
  virtual ~BaseAll() {}

 protected:
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // BASE_ALL_H_

