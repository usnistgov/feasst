
#ifndef FEASST_UTILS_UTILS_IO_H_
#define FEASST_UTILS_UTILS_IO_H_

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <deque>
#include <numeric>
#include <memory>
#include <map>
#include "utils/include/debug.h"
#include "utils/include/utils.h"
#include "utils/include/serialize.h"

using std::cout;
using std::endl;

namespace feasst {

/// Return the feasst install directory as a string.
/// Warning: This does not work with pip or conda Python.
///          Instead, use pyfeasst.forcefield_dir().
inline std::string install_dir() { return std::string(FEASST_DIR_); }

/// Return the feasst version as a string.
inline std::string version() { return std::string(FEASST_VERSION_); }

/// Return the given string with the feasst install directory appended.
/// Warning: This does not work with pip or conda Python.
///          Instead, use pyfeasst.forcefield_dir().
std::string append_install_dir(const char * chars);

/// Return string representation of vector
template<class T>
std::string feasst_str(const std::vector<T> &vec,
  /// use maximum precision
  const bool max_precision = false) {
  std::stringstream ss;
  if (max_precision) ss << MAX_PRECISION;
  for (unsigned int i = 0; i < vec.size(); ++i) {
    ss << vec[i] << " ";
  }
  return ss.str();
}
template<class T>
std::string feasst_str(const std::vector<std::vector<T> > &vec,
  /// use maximum precision
  const bool max_precision = false) {
  std::stringstream ss;
  for (unsigned int i = 0; i < vec.size(); ++i) {
    ss << feasst_str(vec[i], max_precision) << std::endl;
  }
  return ss.str();
}
template<class T>
std::string feasst_str(const std::deque<T> &deq,
  /// use maximum precision
  const bool max_precision = false) {
  std::stringstream ss;
  if (max_precision) ss << MAX_PRECISION;
  for (unsigned int i = 0; i < deq.size(); ++i) {
    ss << deq[i] << " ";
  }
  return ss.str();
}

/// Split a string into a vector of strings via space delimitors.
std::vector<std::string> split(const std::string str);

/// Return phrase with all characters up to the last specialchr removed.
std::string trim(const char* specialchr, const char* phrase,
  /** If 1, remove characters from the left. Otherwise, remove characters
   *  from the right.*/
  int from_left = 1);

/// Same as above except for string phrase.
std::string trim(const char* specialchr, std::string phrase, int from_left = 1);

/// Convert to string with maximum precision
template <typename T>
std::string str(const T a_value) {
  std::ostringstream out;
  out << MAX_PRECISION << a_value;
  return out.str();
}

/// Return the number of spaces in string
int num_spaces(const std::string str);

}  // namespace feasst

#endif  // FEASST_UTILS_UTILS_IO_H_
