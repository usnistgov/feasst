
#ifndef FEASST_UTILS_LET_H_
#define FEASST_UTILS_LET_H_

#include "utils/include/arguments.h"

namespace feasst {

// The actual implementation of this class resides entirely within the arguments
// utilities, but this empty class was created for documentation.
/**
  A Let is a way of defining a string of characters as a variable for reuse on later lines.
 */
class Let {
 public:
  //@{
  /** @name Arguments
    - [var]: a string with any number of spaces, commas, =equal signs, etc.
      The variable name must be enclosed in square brackets, [].
      The value of the variable is then given after the equal sign and can
      contain any number of possible delimitors.
   */
  explicit Let(argtype args = argtype()) {}
  //@}
};

}  // namespace feasst

#endif  // FEASST_UTILS_LET_H_
