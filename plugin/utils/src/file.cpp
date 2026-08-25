#include <sys/stat.h>
#include <sstream>
#include <cstdio>  // std::remove
#include <algorithm>  // count_if
#include "utils/include/file.h"
#include "utils/include/debug.h"

namespace feasst {

void skip_characters(const char comment, std::ifstream &file) {
  std::string line;
  while (file.peek() == comment) {
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

int num_lines(const std::string& file_name) {
  // from https://www.reddit.com/r/cpp_questions/comments/11wlf49/whats_the_most_efficient_way_to_get_the_line/
  std::ifstream file(file_name);
  ASSERT(file.good(), "Could not find file: " << file_name);
  return std::count_if(std::istreambuf_iterator<char>{file}, {}, [](char c) { return c == '\n'; });
}

// from https://stackoverflow.com/a/15119347/31990915
template<typename T1, typename T2>
bool range_equal(T1 first1, T1 last1, T2 first2, T2 last2) {
  while (first1 != last1 && first2 != last2) {
    if (*first1 != *first2) return false;
    ++first1;
    ++first2;
  }
  return (first1 == last1) && (first2 == last2);
}

bool compare_files(const std::string& filename1, const std::string& filename2) {
  std::ifstream file1(filename1);
  std::ifstream file2(filename2);

  std::istreambuf_iterator<char> begin1(file1);
  std::istreambuf_iterator<char> begin2(file2);

  std::istreambuf_iterator<char> end;

  return range_equal(begin1, end, begin2, end);
}

void remove_file(const std::string& filename) {
  if (std::remove(filename.c_str()) != 0) {
    WARN("Error deleting file: " << filename);
  };
}

}  // namespace feasst
