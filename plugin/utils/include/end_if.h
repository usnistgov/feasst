
#ifndef FEASST_UTILS_END_IF_H_
#define FEASST_UTILS_END_IF_H_

#include "utils/include/arguments.h"

namespace feasst {

// The actual implementation of this class resides entirely within the arguments
// utilities, but this empty class was created for documentation.
/**
  EndIf statements terminate the most recent If, and can be nested.
 */
class EndIf {
 public:
  //@{
  /** @name Arguments
   */
  explicit EndIf(argtype args = argtype());
  //@}
};

}  // namespace feasst

#endif  // FEASST_UTILS_END_IF_H_
