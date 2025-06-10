
#ifndef FEASST_UTILS_FOR_H_
#define FEASST_UTILS_FOR_H_

#include "utils/include/arguments.h"

namespace feasst {

// The actual implementation of this class resides entirely within the arguments
// utilities, but this empty class was created for documentation.
/**
  A For loop to repeat the following lines, once for each given variable value,
  until EndFor is reached.
  The variable values are replaced during copying.
  A For loop cannot exist inside of another For loop.
  A For loop must always have a corresponding EndFor.
 */
class For {
 public:
  //@{
  /** @name Arguments
    - [var1]:...:[varN]:val1_0:...:valN_0,...,val1_M:...:valN_M .
      In this example, each line is copied M times while replacing N variables.
      The variables are given by [var1] to [varN] in colon-separated values.
      The values are given by N colon-separated values for each of M copies,
      where copies are comma-separated.
      A For loop requires atleast one variable and value, and ends with EndFor.
   */
  explicit For(argtype args = argtype()) {}
  //@}
};

}  // namespace feasst

#endif  // FEASST_UTILS_FOR_H_
