
#ifndef FEASST_CORE_UTILS_IO_H_
#define FEASST_CORE_UTILS_IO_H_

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <numeric>

using std::cout;
using std::endl;

namespace feasst {

/// Return string representation of vector
template<class T>
std::string feasst_str(const std::vector<T> &vec) {
  std::stringstream ss;
  for (unsigned int i = 0; i < vec.size(); ++i) {
    ss << vec[i] << " ";
  }
  return ss.str();
}
template<class T>
std::string feasst_str(const std::vector<std::vector<T> > &vec) {
  std::stringstream ss;
  for (unsigned int i = 0; i < vec.size(); ++i) {
    ss << feasst_str(vec[i]) << std::endl;
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

}  // namespace feasst

#endif  // FEASST_CORE_UTILS_IO_H_
