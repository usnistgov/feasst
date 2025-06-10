
#ifndef FEASST_UTILS_END_FOR_H_
#define FEASST_UTILS_END_FOR_H_

#include "utils/include/arguments.h"

namespace feasst {

// The actual implementation of this class resides entirely within the arguments
// utilities, but this empty class was created for documentation.
/**
  End a For loop.
 */
class EndFor {
 public:
  //@{
  /** @name Arguments
   */
  explicit EndFor(argtype args = argtype());
  //@}
};

}  // namespace feasst

#endif  // FEASST_UTILS_END_FOR_H_
