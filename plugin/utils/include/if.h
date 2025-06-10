
#ifndef FEASST_UTILS_IF_H_
#define FEASST_UTILS_IF_H_

#include "utils/include/arguments.h"

namespace feasst {

// The actual implementation of this class resides entirely within the arguments
// utilities, but this empty class was created for documentation.
/**
  If statements can detect if the given value is defined or undefined.
  Undefined means that the value is empty (e.g., there is nothing provided).
  Defined means the value is not empty.
 */
class If {
 public:
  //@{
  /** @name Arguments
    - defined=?(value): where the following lines before EndIf or Else are
      used if the given (value) is not empty.
      Otherwise, use any lines between the Else and EndIf.
    - undefined=?(value): same as above, except the following lines are used
      if the given (value) is empty.
   */
  explicit If(argtype args = argtype()) {}
  //@}
};

}  // namespace feasst

#endif  // FEASST_UTILS_IF_H_
