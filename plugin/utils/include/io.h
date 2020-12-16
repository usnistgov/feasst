
#ifndef FEASST_UTILS_IO_H_
#define FEASST_UTILS_IO_H_

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <deque>
#include <iomanip>  // setprecision
#include <limits>  // numeric_limits
#include <map>

namespace feasst {

/// Used to output maximum precision to screen
/// [e.g., INFO(MAX_PRECISION << "energy: " << energy)].
#define MAX_PRECISION std::setprecision(std::numeric_limits<double>::digits10+2)

/// Return the feasst install directory as a string.
/// Warning: This does not work with pip or conda Python.
///          Instead, use pyfeasst.forcefield_dir().
std::string install_dir();

/// Return the feasst version as a string.
std::string version();

/// Return string representation of vector
template<class T>
std::string feasst_str(const std::vector<T> &vec,
  /// use maximum precision
  const bool max_precision = false) {
  std::stringstream ss;
  if (max_precision) ss << MAX_PRECISION;
  for (unsigned int i = 0; i < vec.size(); ++i) {
    ss << vec[i] << ",";
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
    ss << deq[i] << ",";
  }
  return ss.str();
}

/// Print a std::map of a pair of strings in human readable format.
std::string feasst_str(std::map<std::string, std::string> mp);

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

/// Return true if substring is found inside of string.
bool is_found_in(const std::string& str, const std::string& substr);

/// Convert a string to an integer.
int str_to_int(const std::string& str);

/// Convert a string to an double.
double str_to_double(const std::string& str);

/// Convert a string to a boolean.
bool str_to_bool(const std::string& str);

}  // namespace feasst

#endif  // FEASST_UTILS_IO_H_
