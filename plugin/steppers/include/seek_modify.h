
#ifndef FEASST_STEPPERS_SEEK_MODIFY_H_
#define FEASST_STEPPERS_SEEK_MODIFY_H_

#include <string>
#include <vector>
#include "monte_carlo/include/modify.h"

namespace feasst {

class MonteCarlo;

/**
  Find Modify with class name.
 */
class SeekModify {
 public:
  /**
    Return the indices, where the first is mc.modify index.
    If inside ModifyFactory, the second is the index of the factory.
    Otherwise, the second index is -1.
    Only the first match is returned.
   */
  std::vector<int> index(const std::string class_name,
                         const MonteCarlo& mc) const;
};

}  // namespace feasst

#endif  // FEASST_STEPPERS_SEEK_MODIFY_H_
