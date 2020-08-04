#include <sys/stat.h>
#include "utils/include/file.h"
#include "utils/include/debug.h"

namespace feasst {

void skip_characters(const char comment, std::ifstream &file) {
  std::string line;
  std::getline(file, line);
  while (line[0] == comment) {
    std::getline(file, line);
  }
}

bool find(const char* search_string, std::ifstream &file) {
  std::string line;
  std::getline(file, line);
  const int nMax = 1e5;
  int i = 0;
  while ((line.compare(search_string) != 0) &&
         (i != nMax) && (!file.fail())) {
    std::getline(file, line);
    ++i;
    // std::cout << line << std::endl;
  }
  // check if not found
  if (i == nMax || file.fail()) {
    return false;
  }
  return true;
}

bool find(const std::string search_string, std::ifstream &file) {
  return find(search_string.c_str(), file);
}

void find_or_fail(const char* search_string, std::ifstream &file) {
  ASSERT(find(search_string, file), "could not find " << search_string
    << " in file");
}

bool file_exists(const std::string& file_name) {
  struct stat buf;
  if (stat(file_name.c_str(), &buf) != -1) return true;
  return false;
}

void file_backup(const std::string& file_name,
    const std::string append) {
  if (file_exists(file_name)) {
    std::ostringstream ss;
    ss << file_name << append;
    rename(file_name.c_str(), ss.str().c_str());
  }
}

}  // namespace feasst
