
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

/// Return true if file exists.
bool file_exists(const std::string& file_name);

/// If file_name exists, rename file with appended string.
void file_backup(const std::string& file_name,
  const std::string append = ".bak");

/// Return the number of lines in the file.
int num_lines(const std::string& file_name);

/// Return true if the two files are the same
bool compare_files(const std::string& filename1, const std::string& filename2);

/// Remove file
void remove_file(const std::string& filename);

}  // namespace feasst

#endif  // FEASST_UTILS_FILE_H_
