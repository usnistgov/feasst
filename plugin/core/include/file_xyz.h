
#ifndef FEASST_CORE_FILE_XYZ_H_
#define FEASST_CORE_FILE_XYZ_H_

#include "core/include/configuration.h"

namespace feasst {

/// HWH Add html link to XYZ format
class FileXYZ {
 public:
  /// Load the xyz file with file_name into the configuration.
  void load(const char* file_name, Configuration * config) const;

  /// Write the configuration to file_name in xyz format.
  void write(const char* file_name,
             const Configuration& config,
             /// optionally, include group index of configuration to print only
             /// that group.
             int group_index = 0) const;
};

}  // namespace feasst

#endif  // FEASST_CORE_FILE_XYZ_H_
