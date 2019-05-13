
#include <string>
#include <sstream>
#include "utils/include/debug.h"

namespace feasst {

std::string feasst_dir_trim_(const char* file_name) {
  std::stringstream ss;
  ss << FEASST_DIR_;
  const int num_chars = ss.str().size();
  ss.str(file_name);
  std::string file_name_str = ss.str();
  file_name_str.erase(file_name_str.begin(),
                      file_name_str.begin() + num_chars);
  return file_name_str;
}

}  // namespace feasst

