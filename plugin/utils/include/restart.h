
#ifndef FEASST_UTILS_RESTART_H_
#define FEASST_UTILS_RESTART_H_

#include "utils/include/arguments.h"

namespace feasst {

// The actual implementation of this class resides entirely within feasst.cpp
// but this empty class was created for documentation.
/**
  Restart uses a file output by Checkpoint to restart a simulation.
  Checkpoint writes the file every given num_hours, or the file can be written
  immediately by WriteCheckpoint.
 */
class Restart {
 public:
  //@{
  /** @name Arguments
    - checkpoint_file: The file created by Checkpoint or WriteCheckpoint
    - clear_previous_arguments: If true, clear the previous arguments given
      to MonteCarlo, allowing the parsing of new arguments (default: false).
   */
  explicit Restart(argtype args = argtype()) {}
  //@}
};

}  // namespace feasst

#endif  // FEASST_UTILS_RESTART_H_
