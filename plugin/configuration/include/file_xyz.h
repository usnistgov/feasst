
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
      const std::string traj_file_name) {
    std::ofstream vmdf(file_name);
    vmdf << "display projection Orthographic" << endl
      << "color Display Background white" << endl
      << "axes location Off" << endl;
    vmdf << "topo readvarxyz " << trim("/", traj_file_name) << endl;
    vmdf << "mol modstyle 0 0 VDW 1.0000000 120.000000" << endl;
    for (int stype = 0; stype < config.num_site_types(); ++stype) {
      vmdf << "set sel [atomselect top \"name " << stype << "\"]" << endl;
      const double radius = 0.5*config.model_params().sigma().value(stype);
      vmdf << "$sel set radius " << radius << endl;
    }
  }

  void serialize(std::ostream& ostr) const {
    feasst_serialize_version(1, ostr);
  }

  explicit FileVMD(std::istream& istr) {
    feasst_deserialize_version(istr);
  }
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

  void serialize(std::ostream& ostr) const {
    feasst_serialize_version(1, ostr);
    feasst_serialize(group_index_, ostr);
    feasst_serialize(append_, ostr);
  }

  explicit FileXYZ(std::istream& istr) {
    feasst_deserialize_version(istr);
    feasst_deserialize(&group_index_, istr);
    feasst_deserialize(&append_, istr);
  }

 private:
  int group_index_;
  int append_;
};

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_FILE_XYZ_H_
