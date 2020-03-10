
#ifndef FEASST_UTILS_FILE_H_
#define FEASST_UTILS_FILE_H_

#include <string>
#include <fstream>

namespace feasst {

/// Skip all lines beginning with character in file.
void skip_characters(const char comment, std::ifstream &file);

/// Skip all lines in file until reaching a certain line with search_string.
/// Return 1 if found.
bool find(const char* search_string, std::ifstream &file);
bool find(const std::string search_string, std::ifstream &file);

/// Same as find except terminate if search_string is not found.
void find_or_fail(const char* search_string, std::ifstream &file);

}  // namespace feasst

#endif  // FEASST_UTILS_FILE_H_
