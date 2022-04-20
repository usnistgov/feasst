#include <cmath>
#include <iterator>
#include <algorithm>
#include "utils/include/io.h"
#include "utils/include/debug.h"

namespace feasst {

std::string install_dir() {
  return std::string(FEASST_DIR_);
}

std::string version() {
  return std::string(FEASST_VERSION_);
}

// thanks to http://www.cplusplus.com/forum/beginner/87238/
std::vector<std::string> split(const std::string str) {
  std::istringstream buffer(str);
  std::istream_iterator<std::string> begin(buffer), end;
  std::vector<std::string> tokens(begin, end);
  return tokens;
}

std::string trim(const char* specialchr, const char* phrase, int from_left) {
  std::string fs(phrase);

  // find position of last spchr
  std::string spchr(specialchr);
  std::size_t found = fs.find(spchr);
  std::size_t foundprev = 0;
  while (found != std::string::npos) {
    foundprev = found;
    found = fs.find(spchr, found + 1);
  }

  // erase all characters up to spchr
  if (static_cast<int>(foundprev) != 0) {
    if (from_left == 1) {
      fs.erase(fs.begin(), fs.begin() + foundprev + spchr.size());
    } else {
      fs.erase(fs.begin() + foundprev + spchr.size(), fs.end());
    }
  }

  return fs;
}
std::string trim(const char* specialchr, std::string phrase, int from_left) {
  return trim(specialchr, phrase.c_str(), from_left);
}

int num_spaces(const std::string str) {
  return std::count(str.begin(), str.end(), ' ');
}

std::string str(const std::map<std::string, std::string>& mp) {
  std::stringstream ss;
  for (auto const& pair : mp) {
    ss << pair.first << " " << pair.second << " ";
  }
  return ss.str();
}

bool is_found_in(const std::string& str, const std::string& substr) {
  if (str.find(substr) != std::string::npos) return true;
  return false;
}

int str_to_int(const std::string& str) {
  std::stringstream errmsg;
  int intVal = -1;
  errmsg << str << " was " << "expected to be an integer.";
  if (str.find('e') != std::string::npos) {
    const double dble = str_to_double(str);
    const int intt = static_cast<int>(dble);
    ASSERT(std::abs(dble - static_cast<double>(intt)) < 1e-14,
      errmsg.str());
    return intt;
  }
  try {
    intVal = stoi(str);
  } catch (...) {
    FATAL(errmsg.str());
  }
  const double dble = stod(str);
  ASSERT(std::abs(dble - static_cast<double>(intVal)) < 1e-14,
    errmsg.str());
  return intVal;
}

bool str_to_bool(const std::string& str) {
  std::stringstream errmsg;
  errmsg << str << " was expected to be a boolean";
  if (str == "true" || str == "True") {
    return true;
  } else if (str == "false" || str == "False") {
    return false;
  } else if (str == "1") {
    return true;
  } else if (str == "0") {
    return false;
  } else {
    ASSERT(0, errmsg.str());
  }
  return -1;
}

double str_to_double(const std::string& str) {
  DEBUG("str " << str);
  std::stringstream ss(str);
  double double_value = -1;
  ss >> double_value;
  DEBUG("double_value " << MAX_PRECISION << double_value);
  if (ss.fail()) {
    FATAL("given(" << str <<
          ") but was expected to be a double precision number.");
  }
  return double_value;
}

bool replace(const std::string& from, const std::string& to, std::string * str) {
  size_t start_pos = str->find(from);
  if (start_pos == std::string::npos) {
    return false;
  }
  str->replace(start_pos, from.length(), to);
  return true;
}

}  // namespace feasst
