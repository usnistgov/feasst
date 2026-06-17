#include <string>
#include <iostream>
#ifdef _OPENMP
  #include <omp.h>
#endif  // _OPENMP
#include "utils/include/io.h"
#include "utils/include/debug.h"
#include "utils/include/definitions.h"

namespace feasst {

std::string feasst_dir_trim_(const char* file_name) {
  //std::cout << "filename " << file_name << std::endl;
  std::string text = file_name;
  size_t pos = text.find("/feasst/");
  if (pos != std::string::npos) {
    text.erase(0, pos); // Erase from index 0 up to pos
  }
  return text;
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
