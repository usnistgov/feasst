
#ifndef FEASST_CONFIGURATION_FILE_VMD_H_
#define FEASST_CONFIGURATION_FILE_VMD_H_

#include <map>
#include <string>
#include <fstream>
#include "configuration/include/configuration.h"
#include "configuration/include/visit_configuration.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  Print a vmd script to view an xyz file via Bash: "vmd -e file.vmd"

  For macOS, consider using an alias in your .bash_profile such as:
  alias vmd='/Applications/VMD\ 1.9.4a51-x86_64-Rev9.app/Contents/MacOS/startup.command'
 */
class FileVMD {
 public:
  //@{
  /** @name Arguments
   */

  /**
    args:
    - min_sigma: minimum sigma for printing radius in VMD (default: 0.1).
   */
  explicit FileVMD(argtype args = argtype());
  explicit FileVMD(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  // Place holder for derived classes (e.g. FileVMDPatch)
  virtual void get_params(const Configuration& config,
    const int site_type,
    double * radius,
    double * distance,
    int * center_index) const;

  void write(const std::string file_name,
    const Configuration& config,
    const std::string traj_file_name) const;

  void serialize(std::ostream& ostr) const;
  explicit FileVMD(std::istream& istr);
  virtual ~FileVMD() {}

  //@}
 private:
  double min_sigma_;
};

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_FILE_VMD_H_
