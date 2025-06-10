
#ifndef FEASST_UTILS_ELSE_H_
#define FEASST_UTILS_ELSE_H_

#include "utils/include/arguments.h"

namespace feasst {

// The actual implementation of this class resides entirely within the arguments
// utilities, but this empty class was created for documentation.
/**
  Else statements use the following lines until EndIf when the preceeding If is
  not true.
 */
class Else {
 public:
  //@{
  /** @name Arguments
   */
  explicit Else(argtype args = argtype());
  //@}
};

}  // namespace feasst

#endif  // FEASST_UTILS_ELSE_H_
