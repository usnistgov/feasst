
#ifndef FEASST_NETCDF_FILE_NETCDF_H_
#define FEASST_NETCDF_FILE_NETCDF_H_

#include <memory>
#include <vector>
#include <string>
#include <fstream>
#include <netcdf>
#include "configuration/include/configuration.h"
#include "configuration/include/visit_configuration.h"

namespace feasst {

/**
 */
class FileNETCDF {
 public:
  // HWH: These are all the same arguments as FileXYZ.
  // HWH: Perhaps make a File base class to share them.
  // HWH: Always append
  /**
    args:
    - group_index: print the coordinates of this group index only (default: 0).
    - group: name of group defined within system (default: "").
    - euler: if true, print Euler angles (default: "false").
    - float_precision: if true, use floating point precision for coordinates.
      Otherwise, use double precision (default: false).
   */
  explicit FileNETCDF(argtype args = argtype());
  explicit FileNETCDF(argtype * args);

  /// Initialize the dimensions.
  void initialize(const Configuration& config);

  /// Write the configuration.
  /// If the simulation is 2D, simply write z as 0.
  void write(const Configuration& config);

  /**
    Load the xyz file with file_name into the configuration.
    Note that this function does not read the domain tilts xy, xz and yz.
    This function also does not read the site types.
    Thus, particles should be added to the system in the desired order.
    If no particles in the configuration, use the first particle type.
   */
  void load(const std::string file_name, Configuration * config) const;

  /// As above, but can load frame by frame.
  void load_frame(std::ifstream& xyz, Configuration * config) const;

  void serialize(std::ostream& ostr) const;
  explicit FileNETCDF(std::istream& istr);

 private:
  int group_index_;
  std::string group_;
  bool euler_;
  bool float_precision_;
  int num_config_;

  const std::shared_ptr<netCDF::NcFile> file_;
  std::vector<netCDF::NcDim> xyz_dims_;
  std::vector<netCDF::NcDim> domain_dims_;
  std::vector<netCDF::NcDim> type_dims_;
  netCDF::NcDim record_dim_;
  netCDF::NcVar nc_domain_;
  netCDF::NcVar nc_tilt_xy_;
  netCDF::NcVar nc_tilt_xz_;
  netCDF::NcVar nc_tilt_yz_;
  netCDF::NcVar nc_site_types_;
  netCDF::NcVar nc_num_sites_;
  netCDF::NcVar nc_xyzs_;

  int gindex_(const Configuration& config) const {
    int gindex = group_index_;
    if (!group_.empty()) {
      gindex = config.group_index(group_);
    }
    return gindex;
  }

  const std::shared_ptr<netCDF::NcFile> file_parse_(argtype * args) {
    const std::string file_name = str("file_name", args, "");
    if (!file_name.empty()) {
      return std::make_shared<netCDF::NcFile>(file_name,
                                              netCDF::NcFile::replace);
    }
    return NULL;
  }
};

// Utility class to print XYZ files from selection.
class PrinterNETCDF : public LoopConfigOneBody {
 public:
  PrinterNETCDF(
             const bool euler,
             const bool float_precision,
             std::vector<int> * site_types,
             double * xyzs,
             float * xyzs_f);
  void work(const Site& site,
    const Configuration& config,
    const LoopDescriptor& data) override;
 private:
  bool euler_;
  bool float_precision_;
  std::vector<int> * site_types_;
  double * xyzs_;
  float * xyzs_f_;
  int site_index_ = 0;
};

}  // namespace feasst

#endif  // FEASST_NETCDF_FILE_NETCDF_H_
