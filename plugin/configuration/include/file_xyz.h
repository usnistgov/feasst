
#ifndef FEASST_CONFIGURATION_FILE_XYZ_H_
#define FEASST_CONFIGURATION_FILE_XYZ_H_

#include <fstream>
#include <map>
#include <memory>
#include <string>
#include "configuration/include/configuration.h"
#include "configuration/include/visit_configuration.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

// Utility class to print XYZ files from selection.
class PrinterXYZ : public LoopConfigOneBody {
 public:
  PrinterXYZ(std::shared_ptr<std::ofstream> file,
             const bool euler,
             const int num_places = 8);
  void work(const Site& site,
    const Configuration& config,
    const LoopDescriptor& data) override;
 private:
  int num_places_;
  bool euler_;
  std::shared_ptr<std::ofstream> file_;
};

// Note HWH: best not to assume this second line format by default
/**
  The XYZ file format has no formal standard: https://en.wikipedia.org/wiki/XYZ_file_format
  It is important to remember that FEASST has its own variant.

  The first line is the number of sites, n.
  The second line is of the format [id lx ly lz xy xz yz], where id is a
  placeholder for order parameters or macrostates,
  lx is the box length in the x direction,
  ly is the box length in the y direction,
  lz is the box length in the z direction,
  xy is the xy domain tilt (see Domain),
  xz is the xz domain tilt, and yz is the yz domain tilt.
  If lz == 0, then the domain is set to two dimensions.
  The following n lines are in the format [id x y z],
  where id is the unique site type and x, y, z are the Cartesian coordinates.

  Although FileXYZ writes the domain tilt factors, FileXYZ does not read them.
 */
class FileXYZ {
 public:
  //@{
  /** @name Arguments
   */

  /**
    args:
    - group_index: print the coordinates of this group index only (default: 0).
    - group: name of group defined within system (default: "").
    - append: append file output if set to true.
      Do not append if false (default: "false").
    - euler: if true, print Euler angles (default: "false").
   */
  explicit FileXYZ(argtype args = argtype());
  explicit FileXYZ(argtype * args);

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

  /// As above, but can load frame by frame.
  /// Return false if the file as at the end.
  bool load_frame(std::ifstream& xyz, Configuration * config) const;

  /// Write the configuration to file_name in xyz format.
  /// If the simulation is 2D, simply writes z as 0.
  void write(const std::string file_name,
             const Configuration& config,
             /// Number of decimal places
             const int num_decimal_places = 8) const;

  void write_for_vmd(const std::string file_name,
             const Configuration& config) const;

  void serialize(std::ostream& ostr) const;
  explicit FileXYZ(std::istream& istr);
  ~FileXYZ();

  //@}
 private:
  int group_index_;
  std::string group_;
  bool append_, euler_;
};

inline std::shared_ptr<FileXYZ> MakeFileXYZ(argtype args = argtype()) {
  return std::make_shared<FileXYZ>(args);
}

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_FILE_XYZ_H_
