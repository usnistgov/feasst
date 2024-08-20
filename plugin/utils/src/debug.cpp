#include <string>
#include <iostream>
#ifdef _OPENMP
  #include <omp.h>
#endif  // _OPENMP
#include "utils/include/io.h"
#include "utils/include/debug.h"

namespace feasst {

std::string feasst_dir_trim_(const char* file_name) {
  std::stringstream ss;
  ss << FEASST_DIR_ "/";
  const int num_chars = ss.str().size();
  ss.str(file_name);
  std::string file_name_str = ss.str();
  file_name_str.erase(file_name_str.begin(),
                      file_name_str.begin() + num_chars);
  return file_name_str;
}

void feasst_macro_output(const std::string& name, std::string message) {
#ifdef _OPENMP
  message = "#" + name + str(omp_get_thread_num()) + message;
#else
  message = "#" + name + message;
#endif  // _OPENMP
  std::cout << message << std::endl;
}

std::string feasst_omp_thread() {
# ifdef _OPENMP
  return str(omp_get_thread_num());
# else
  return "";
# endif  // _OPENMP
}

}  // namespace feasst
