
#ifndef FEASST_PATCH_FILE_XYZ_PATCH_H_
#define FEASST_PATCH_FILE_XYZ_PATCH_H_

#include <string>
#include <fstream>
#include <sstream>
#include "utils/include/arguments.h"
#include "configuration/include/configuration.h"
#include "configuration/include/visit_configuration.h"
#include "configuration/include/file_xyz.h"

namespace feasst {

/**
  Print a vmd script to view an xyz file via command "vmd -e file.vmd"

  Make a patch by inscribing sphere of radius r2 inside bead of radius r1=sig.
  The distance between center of sphere and bead is given by "d."
  For given patch angle, use law of cosines to derive r2 via
  constraint r1 + eps = r2 + d, where eps is small, to avoid clipping.
 */
class FileVMDPatch : public FileVMD {
 public:
  FileVMDPatch() : FileVMD() {}

  void get_params(const Configuration& config,
                  const int site_type,
                  double * radius,
                  double * distance,
                  int * center_index) const override;

  void serialize(std::ostream& ostr) const;
  explicit FileVMDPatch(std::istream& istr);
};

// Utility class to print XYZ files from selection.
class PrinterXYZPatch : public LoopConfigOneBody {
 public:
  PrinterXYZPatch(std::shared_ptr<std::ofstream> file, const int num_places = 8);
  void work(const Site& site,
      const Configuration& config,
      const LoopDescriptor& data) override;

 private:
  int num_places_;
  std::shared_ptr<std::ofstream> file_;
  FileVMDPatch vmd_;
};

/**
  Similar to FileXYZ, except for patches.
 */
class FileXYZPatch {
 public:
  /**
    args:
    - group_index: print the coordinates of this group index only (default: 0).
    - append: append file output if set to true.
      Do not append if false (default: "false").
   */
  explicit FileXYZPatch(argtype args = argtype());
  explicit FileXYZPatch(argtype * args);

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
  explicit FileXYZPatch(std::istream& istr);

 private:
  int group_index_;
  bool append_;
};

inline std::shared_ptr<FileXYZPatch> MakeFileXYZPatch(argtype args = argtype()) {
  return std::make_shared<FileXYZPatch>(args);
}

}  // namespace feasst

#endif  // FEASST_PATCH_FILE_XYZ_PATCH_H_
