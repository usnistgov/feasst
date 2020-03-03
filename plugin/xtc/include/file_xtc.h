
#ifndef FEASST_XTC_FILE_XTC_H_
#define FEASST_XTC_FILE_XTC_H_

#include <string>
#include "xdrfile.h"
#include "xdrfile_xtc.h"
#include "configuration/include/configuration.h"

namespace feasst {

/**
  XTC files are binary compressed coordinate and box trajectories:
  http://www.gromacs.org/Developer_Zone/Programming_Guide/XTC_Library

  This implementation assumes 3 dimensions.
  This implementation also currently assumes a cuboid domain.
  HWH add triclinic.

  To read an XDRFILE in C++, you must open it as follows:

  XDRFILE * xdr_file = xdrfile_open(file_name_c_str, "r");

  To write an XDRFILE in C++, you must open it as follows:

  XDRFILE * xdr_file = xdrfile_open(file_name_c_str, "w");
 */
class FileXTC {
 public:
  /// Load XTC file into configuration. Return 0 if successful.
  int load(XDRFILE * file, Configuration * config,
    /// provide the file_name to check number of atoms. Skip if empty.
    std::string file_name = "");

  /// Write the configuration to XTC file
  void write(XDRFILE * file, const Configuration& config);
};

}  // namespace feasst

#endif  // FEASST_XTC_FILE_XTC_H_
