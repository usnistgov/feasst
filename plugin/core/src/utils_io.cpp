
#include <iterator>
#include <algorithm>
#include "core/include/utils_io.h"

namespace feasst {

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
  if (int(foundprev) != 0) {
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

void feasst_serialize(const std::string str, std::ostream& ostr) {
  ASSERT(num_spaces(str) == 0, "no spaces in serialized string");
  if (str.empty()) {
    ostr << "0 ";
  } else {
    ostr << "1 ";
    ostr << str << " ";
  }
}

void feasst_deserialize(std::string * str, std::istream& istr) {
  int empty;
  istr >> empty;
  if (empty != 0) {
    istr >> *str;
  }
}

void feasst_serialize(const double val, std::ostream& ostr) {
  ostr << MAX_PRECISION << val << " ";
}

void feasst_deserialize(double * val, std::istream& istr) {
  istr >> *val;
}

void feasst_serialize(const int val, std::ostream& ostr) {
  ostr << val << " ";
}

void feasst_deserialize(int * val, std::istream& istr) {
  istr >> *val;
}

void feasst_serialize(const bool val, std::ostream& ostr) {
  ostr << val << " ";
}

void feasst_deserialize(bool * val, std::istream& istr) {
  int tmp;
  istr >> tmp;
  *val = tmp;
}

void feasst_serialize_version(const int version, std::ostream& ostr) {
  ostr << version << " ";
}


int feasst_deserialize_version(std::istream& istr) {
  int version;
  istr >> version;
  return version;
}

}  // namespace feasst
