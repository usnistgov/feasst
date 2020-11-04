#include <limits>
#include "utils/include/debug.h"
#include "utils/include/serialize.h"

namespace feasst {

void feasst_serialize(const std::string str, std::ostream& ostr) {
  ASSERT(num_spaces(str) == 0, "no spaces allowed in serialized string("
    << str << ")");
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

void feasst_serialize_endcap(const std::string name, std::ostream& ostr) {
  ostr << "End" << name << " ";
}

void feasst_deserialize_endcap(const std::string name, std::istream& istr) {
  std::string read_name;
  istr >> read_name;
  ASSERT("End" + name == read_name, "There is a problem with serialization. "
    << "The endcap was expected to be: " << name << " "
    << "but was found to be: " << read_name);
}

void feasst_deserialize(double * val, std::istream& ostr) {
  std::string valstr;
  ostr >> valstr;
  if (valstr == "inf") {
    *val = 2*std::numeric_limits<double>::max();
  } else {
    try {
      *val = std::stod(valstr);
    } catch (...) {
      *val = 2*std::numeric_limits<double>::max();
    }
  }
}

void feasst_serialize(const long double& val, std::ostream& ostr) {
  ostr << std::setprecision(std::numeric_limits<long double>::digits10+2)
       << val << " ";
}


void feasst_deserialize(long double * val, std::istream& ostr) {
  std::string valstr;
  ostr >> valstr;
  if (valstr == "inf") {
    *val = 2*std::numeric_limits<long double>::max();
  } else {
    try {
      *val = std::stold(valstr);
    } catch (...) {
      *val = 2*std::numeric_limits<long double>::max();
    }
  }
}

void feasst_serialize(const std::vector<double>& vector, std::ostream& ostr) {
  ostr << MAX_PRECISION;
  ostr << vector.size() << " ";
  for (const double& element : vector) {
    feasst_serialize(element, ostr);
  }
}

void feasst_deserialize(std::vector<double> * vector, std::istream& istr) {
  int num;
  istr >> num;
  vector->resize(num);
  for (int index = 0; index < num; ++index) {
    feasst_deserialize(&(*vector)[index], istr);
  }
}

void feasst_serialize(const std::vector<long double>& vector, std::ostream& ostr) {
  ostr << MAX_PRECISION;
  ostr << vector.size() << " ";
  for (const long double& element : vector) {
    feasst_serialize(element, ostr);
  }
}

void feasst_deserialize(std::vector<long double> * vector, std::istream& istr) {
  int num;
  istr >> num;
  vector->resize(num);
  for (int index = 0; index < num; ++index) {
    feasst_deserialize(&(*vector)[index], istr);
  }
}

void feasst_serialize(const argtype& args, std::ostream& ostr) {
  feasst_serialize(args.size(), ostr);
  for (const auto& pair : args) {
    feasst_serialize(pair.first, ostr);
    feasst_serialize(pair.second, ostr);
  }
}

void feasst_deserialize(argtype * args, std::istream& istr) {
  int num;
  feasst_deserialize(&num, istr);
  std::string first, second;
  for (int i = 0; i < num; ++i) {
    feasst_deserialize(&first, istr);
    feasst_deserialize(&second, istr);
    args->insert({first, second});
  }
}

}  // namespace feasst
