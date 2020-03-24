
#ifndef FEASST_CONFIGURATION_FILE_XYZ_H_
#define FEASST_CONFIGURATION_FILE_XYZ_H_

#include <string>
#include <fstream>
#include <sstream>
#include "configuration/include/configuration.h"

namespace feasst {

class FileVMD {
 public:
  FileVMD() {}

  void write(const std::string file_name,
      const Configuration& config,
      const std::string traj_file_name);

  void serialize(std::ostream& ostr) const;
  explicit FileVMD(std::istream& istr);
};

/// HWH Add html link to XYZ format
/// Note that the load() function reads the box length from the second line
/// according to the format [id lx ly lz]
/// Note HWH: best not to assume this second line format by default
class FileXYZ {
 public:
  FileXYZ() {
    set_group_index();
    set_append();
  }

  /// Load the xyz file with file_name into the configuration.
  /// If no particles in the configuration, try to use the first particle type.
  void load(const std::string file_name, Configuration * config) const;

  /// Write the configuration to file_name in xyz format.
  void write(const std::string file_name,
             const Configuration& config) const;

  void write_for_vmd(const std::string file_name,
             const Configuration& config) const;

  void set_group_index(const int index = 0) { group_index_ = index; }

  /// By default, do not append
  void set_append(const int append = 0) { append_ = append; }

  void serialize(std::ostream& ostr) const;
  explicit FileXYZ(std::istream& istr);

 private:
  int group_index_;
  int append_;
};

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_FILE_XYZ_H_
