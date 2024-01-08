
#ifndef FEASST_SPHEROCYLINDER_FILE_XYZ_SPHEROCYLINDER_H_
#define FEASST_SPHEROCYLINDER_FILE_XYZ_SPHEROCYLINDER_H_

#include <string>
#include <fstream>
#include <sstream>
#include "utils/include/arguments.h"
#include "configuration/include/configuration.h"
#include "configuration/include/visit_configuration.h"
#include "configuration/include/file_xyz.h"
#include "configuration/include/file_vmd.h"

namespace feasst {

/**
  Print a vmd script to view an xyz file via command "vmd -e file.vmd"

  Make a spherocylinder by placing a sphere on each endcap and then connecting
  them in VMD using the CPK representation.
 */
class FileVMDSpherocylinder : public FileVMD {
 public:
  FileVMDSpherocylinder() : FileVMD() {}

  void get_params(const Configuration& config,
                  const int site_type,
                  double * radius,
                  double * length,
                  int * center_index) const override;

  void serialize(std::ostream& ostr) const;
  explicit FileVMDSpherocylinder(std::istream& istr);
};

// Utility class to print XYZ files from selection.
class PrinterXYZSpherocylinder : public LoopConfigOneBody {
 public:
  PrinterXYZSpherocylinder(std::shared_ptr<std::ofstream> file, const int num_places, const int expand_factor);
  void work(const Site& site,
      const Configuration& config,
      const LoopDescriptor& data) override;

 private:
  int num_places_;
  int expand_factor_;
  std::shared_ptr<std::ofstream> file_;
  FileVMDSpherocylinder vmd_;
};

/**
  Similar to FileXYZ, except for spherocylinders.
 */
class FileXYZSpherocylinder {
 public:
  //@{
  /** @name Arguments
    - expand_factor: represent spherocylinder as this many spheres (default: 10).
    - group_index: print the coordinates of this group index only (default: 0).
    - append: append file output if set to true.
      Do not append if false (default: "false").
   */
  explicit FileXYZSpherocylinder(argtype args = argtype());
  explicit FileXYZSpherocylinder(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /**
    Load the xyz file with file_name into the configuration.
    Note that this function does not read the domain tilts xy, xz and yz.
    This function also does not read the site types.
    Thus, particles should be added to the system in the desired order.
    If no particles in the configuration, use the first particle type.
   */
  void load(const std::string file_name, Configuration * config) const;

  /// Write the configuration to file_name in xyz format.
  /// If the simulation is 2D, simply writes z as 0.
  void write(const std::string file_name,
             const Configuration& config,
             /// Number of decimal places
             const int num_decimal_places = 8) const;

  void write_for_vmd(const std::string file_name,
             const Configuration& config) const;
  bool append() const { return append_; }
  void serialize(std::ostream& ostr) const;
  explicit FileXYZSpherocylinder(std::istream& istr);

  //@}
 private:
  int expand_factor_;
  int group_index_;
  bool append_;
};

inline std::shared_ptr<FileXYZSpherocylinder> MakeFileXYZSpherocylinder(argtype args = argtype()) {
  return std::make_shared<FileXYZSpherocylinder>(args);
}

}  // namespace feasst

#endif  // FEASST_SPHEROCYLINDER_FILE_XYZ_SPHEROCYLINDER_H_
