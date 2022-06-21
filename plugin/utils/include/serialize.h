
#ifndef FEASST_UTILS_SERIALIZE_H_
#define FEASST_UTILS_SERIALIZE_H_

#include <string>
#include <sstream>
#include <vector>
#include <deque>
#include <memory>
#include <map>
#include "utils/include/arguments.h"
#include "utils/include/debug.h"
#include "utils/include/io.h"

namespace feasst {

/**
  Utility functions for serialization of objects into human-readable character
  streams.
  Each serialization function has an accompying deserialization function which
  performs the inverse operation.

  Note that some deserialization functions of derived objects do not work as
  the commented templates.
  But when the source code is copied to the specific object's deserialization
  function, then it does work, for unknown reasons.

  For this function, serialize a boolean value.
 */
void feasst_serialize(const bool val, std::ostream& ostr);

/// Deserialize bool.
void feasst_deserialize(bool * val, std::istream& istr);

/// Serialize string.
void feasst_serialize(const std::string str, std::ostream& ostr);

/// Deserialize string.
void feasst_deserialize(std::string * str, std::istream& istr);

/// Serialize double. Handle zero.
void feasst_serialize(const double val, std::ostream& ostr);

/// Deserialize double. Handle inf.
void feasst_deserialize(double * val, std::istream& istr);

/// Serialize long double.
void feasst_serialize(const long double& val, std::ostream& ostr);

/// Deserialize long double. Handle inf.
void feasst_deserialize(long double * val, std::istream& istr);

/// Serialize argtype.
void feasst_serialize(const argtype& args, std::ostream& ostr);

/// Deserialize argtype.
void feasst_deserialize(argtype * args, std::istream& istr);

/// Serialize generic to full precision.
template <typename T>
void feasst_serialize(const T& val, std::ostream& ostr) {
  ostr << MAX_PRECISION << val << " ";
}

/// Deserialize generic
template <typename T>
void feasst_deserialize(T * val, std::istream& istr) {
  istr >> *val;
}

/// Serialize object version
void feasst_serialize_version(const int version, std::ostream& ostr);

/// Deserialize object version
int feasst_deserialize_version(std::istream& istr);

/// Serialize pair of T1 and T2 vector.
template <typename T1, typename T2>
void feasst_serialize(const std::pair<T1, std::vector<T2> >& container,
  std::ostream& ostr) {
  feasst_serialize(container.first, ostr);
  ostr << container.second.size() << " ";
  for (const T2& element : container.second) {
    feasst_serialize(element, ostr);
  }
}

/// deserialize pair of T1 and T2 vector.
template <typename T1, typename T2>
void feasst_deserialize(std::pair<T1, std::vector<T2> > * container,
  std::istream& istr) {
  feasst_deserialize(&container->first, istr);
  int num;
  istr >> num;
  container->second.resize(num);
  for (int index = 0; index < num; ++index) {
    feasst_deserialize(&(container->second)[index], istr);
  }
}

/// Serialize pair of T1 and T2.
template <typename T1, typename T2>
void feasst_serialize(const std::pair<T1, T2>& container,
  std::ostream& ostr) {
  feasst_serialize(container.first, ostr);
  feasst_serialize(container.second, ostr);
}

/// deserialize pair of T1 and T2.
template <typename T1, typename T2>
void feasst_deserialize(std::pair<T1, T2> * container,
  std::istream& istr) {
  feasst_deserialize(&container->first, istr);
  feasst_deserialize(&container->second, istr);
}

/// Serialize the 1D vector.
template <typename T>
void feasst_serialize(const std::vector<T>& vector, std::ostream& ostr) {
  ostr << vector.size() << " ";
  for (const T& element : vector) {
    feasst_serialize(element, ostr);
  }
}

/// Deserialize the 1D vector.
template <typename T>
void feasst_deserialize(std::vector<T> * vector, std::istream& istr) {
  int num;
  istr >> num;
  vector->resize(num);
  for (int index = 0; index < num; ++index) {
    feasst_deserialize(&(*vector)[index], istr);
  }
}

/// Serialize 1D vector of doubles
void feasst_serialize(const std::vector<double>& vector, std::ostream& ostr);

/// Deserialize 1D vector of doubles.
void feasst_deserialize(std::vector<double> * vector, std::istream& istr);

/// Serialize 1D vector of long doubles
void feasst_serialize(const std::vector<long double>& vector, std::ostream& ostr);

/// Deserialize 1D vector of long doubles.
void feasst_deserialize(std::vector<long double> * vector, std::istream& istr);

/// Serialize the 1D deque.
template <typename T>
void feasst_serialize(const std::deque<T>& deque, std::ostream& ostr) {
  ostr << MAX_PRECISION;
  ostr << deque.size() << " ";
  for (const T& element : deque) {
    ostr << element << " ";
  }
}

/// Deserialize the 1D deque.
template <typename T>
void feasst_deserialize(std::deque<T> * deque, std::istream& istr) {
  int num;
  istr >> num;
  deque->resize(num);
  for (int index = 0; index < num; ++index) {
    istr >> (*deque)[index];
  }
}

/// Deserialize the boolean 1D vector.
inline void feasst_deserialize(std::vector<bool> * vector, std::istream& istr) {
  int num;
  istr >> num;
  vector->resize(num);
  for (int index = 0; index < num; ++index) {
    int tmp;
    istr >> tmp;
    (*vector)[index] = tmp;
  }
}

/// Deserialize the boolean 2D vector.
inline void feasst_deserialize(std::vector<std::vector<bool> > * vector,
                               std::istream& istr) {
  int dim1;
  istr >> dim1;
  vector->resize(dim1);
  for (int index1 = 0; index1 < dim1; ++index1) {
    int dim2;
    istr >> dim2;
    (*vector)[index1].resize(dim2);
    for (int index2 = 0; index2 < dim2; ++index2) {
      int tmp;
      istr >> tmp;
      (*vector)[index1][index2] = tmp;
    }
  }
}

/// Serialize feasst object
template <typename T>
void feasst_serialize_fstobj(const T& obj, std::ostream& ostr) {
  obj.serialize(ostr);
}

/// Deserialize feasst object
template <typename T>
void feasst_deserialize_fstobj(T * obj, std::istream& istr) {
  *obj = T(istr);
}

/// Serialize vector of feasst objects
template <typename T>
void feasst_serialize_fstobj(const std::vector<T>& vector, std::ostream& ostr) {
  ostr << vector.size() << " ";
  for (const T& element : vector) {
    element.serialize(ostr);
  }
}

/// Deserialize vector of feasst objects
template <typename T>
void feasst_deserialize_fstobj(std::vector<T> * vector, std::istream& istr) {
  int dim1;
  istr >> dim1;
  // vector->resize(dim1);  // push_back avoids constructor with required args.
  for (int index = 0; index < dim1; ++index) {
    vector->push_back(T(istr));
  }
}

/// Serialize 2D vector of feasst objects
template <typename T>
void feasst_serialize_fstobj(const std::vector<std::vector<T> >& vector,
    std::ostream& ostr) {
  ostr << vector.size() << " ";
  for (const std::vector<T>& element : vector) {
    feasst_serialize_fstobj(element, ostr);
  }
}

/// Deserialize 2D vector of feasst objects
template <typename T>
void feasst_deserialize_fstobj(std::vector<std::vector<T> > * vector,
    std::istream& istr) {
  int dim1;
  istr >> dim1;
  for (int index = 0; index < dim1; ++index) {
    std::vector<T> vec;
    feasst_deserialize_fstobj(&vec, istr);
    vector->push_back(vec);
  }
}

/// Serialize double stored as shared pointer
inline void feasst_serialize_sp(const std::shared_ptr<double> ptr,
    std::ostream& ostr) {
  if (ptr) {
    ostr << "1 ";
    feasst_serialize(*ptr, ostr);
  } else {
    ostr << "0 ";
  }
}

// HWH for unknown reasons, this function template does not work.
// /// Deserialize double stored as shared pointer
// inline void feasst_deserialize(std::shared_ptr<double> ptr,
//     std::istream& istr) {
//   int existing;
//   istr >> existing;
//   if (existing != 0) {
//     double value;
//     istr >> value;
//     ptr = std::make_shared<double>(value);
//   }
// }

/// Serialize int stored as shared pointer
inline void feasst_serialize_sp(const std::shared_ptr<int> ptr,
    std::ostream& ostr) {
  if (ptr) {
    ostr << "1 ";
    feasst_serialize(*ptr, ostr);
  } else {
    ostr << "0 ";
  }
}

// HWH for unknown reasons, this function template does not work.
// /// Deserialize int stored as shared pointer
// inline void feasst_deserialize(std::shared_ptr<int> ptr,
//     std::istream& istr) {
//   int existing;
//   istr >> existing;
//   if (existing != 0) {
//     int value;
//     istr >> value;
//     ptr = std::make_shared<int>(value);
//   }
// }

/// Serialize feasst object stored as shared pointer
template <typename T>
void feasst_serialize(const std::shared_ptr<T> ptr, std::ostream& ostr) {
  if (ptr) {
    ostr << "1 ";
    ptr->serialize(ostr);
  } else {
    ostr << "0 ";
  }
}

/// Deserialize feasst object stored as shared pointer
template <typename T>
void feasst_deserialize(std::shared_ptr<T> ptr, std::istream& istr) {
  int existing;
  istr >> existing;
  if (existing != 0) {
    ptr = std::make_shared<T>(istr);
  }
}

/// Serialize feasst derived object stored as shared pointer
template <typename T>
void feasst_serialize_fstdr(std::shared_ptr<T> ptr, std::ostream& ostr) {
  feasst_serialize(ptr, ostr);
}

// HWH for unknown reasons, this function template does not work.
// /// Deserialize feasst derived object stored as shared pointer
// template <typename T>
// void feasst_deserialize_fstdr(std::shared_ptr<T> ptr, std::istream& istr) {
//   int existing;
//   istr >> existing;
//   if (existing != 0) {
//     ptr = ptr->deserialize(istr);
//   }
// }

/// Serialize vector of shared pointers of feasst objects
template <typename T>
void feasst_serialize(const std::vector<std::shared_ptr<T> >& vector,
    std::ostream& ostr) {
  ostr << vector.size() << " ";
  for (const std::shared_ptr<T> element : vector) {
    feasst_serialize(element, ostr);
  }
}

// HWH this didn't work in Accumulator. Depreciate it?
/// Deserialize vector of shared pointers of feasst objects
template <typename T>
void feasst_deserialize(std::vector<std::shared_ptr<T> > * vector,
    std::istream& istr) {
  int dim1;
  istr >> dim1;
  vector->resize(dim1);
  for (int index = 0; index < dim1; ++index) {
    feasst_deserialize((*vector)[index], istr);
  }
}

/// Serialize vector of shared pointers of feasst derived objects
template <typename T>
void feasst_serialize_fstdr(const std::vector<std::shared_ptr<T> >& vector,
    std::ostream& ostr) {
  ostr << vector.size() << " ";
  for (std::shared_ptr<T> element : vector) {
    feasst_serialize_fstdr(element, ostr);
  }
}

//  HWH for unknown reasons, this function template does not work.
// /// Deserialize vector of shared pointers of feasst derived objects
// template <typename T>
// void feasst_deserialize_fstdr(std::vector<std::shared_ptr<T> > * vector,
//     std::istream& istr) {
//   int dim1;
//   istr >> dim1;
//   vector->resize(dim1);
//   for (int index = 0; index < dim1; ++index) {
//     feasst_deserialize_fstdr((*vector)[index], istr);
//   }
// }

/// Return a shared pointer to the base class of model after construction of
/// the full derived class.
/// see https://isocpp.org/wiki/faq/serialization
template <typename T>
std::shared_ptr<T> template_deserialize(
    std::map<std::string, std::shared_ptr<T> >& map,
    std::istream& istr,
    /// Rewind istr position to read class name again (default: false).
    bool rewind = false) {
  std::string class_name;
  int pos = istr.tellg();  // record position before reading
  istr >> class_name;      // read class name

  // rewind position so constructors can reread class name.
  if (rewind) {
    istr.seekg(pos, istr.beg);  // rewind to before reading the class name.
  }
  DEBUG("deserializing: " << class_name << " rewind? " << rewind);
  if (map.count(class_name) == 0) {
    FATAL("The class name \"" << class_name << "\" "
    << "is not recognized during deserialization. "
    << "If the above class name is empty, there was a mis-match in stream. "
    << "Perhaps the plugin was not included during compilation. "
    << "If that's not it, its likely due to the lack of a static mapper "
    << "which is typically implemented within the cpp file. "
    << "In rare cases, the absence of a constructor implementation inside "
    << "the cpp file possibly leads optimization to ignore the mapper.");
  }
  std::shared_ptr<T> obj = map[class_name]->create(istr);
  DEBUG("obj " << obj);
  return obj;
}

/// Return a deep copy of a feasst derived class object.
/// This is implemented via serialization/deserialization.
template <typename T>
std::shared_ptr<T> deep_copy_derived(std::shared_ptr<T> object) {
  std::stringstream ss;
  object->serialize(ss);
  return object->deserialize(ss);
}
template <typename T>
std::shared_ptr<T> deep_copy_derived(T * object) {
  std::stringstream ss;
  object->serialize(ss);
  return object->deserialize(ss);
}

/// Return a deep copy.
/// This is implemented via serialization/deserialization.
template <typename T>
T deep_copy(const T& object) {
  std::stringstream ss;
  object.serialize(ss);
  return T(ss);
}

/// End class serialization with this notification to aid debugging.
void feasst_serialize_endcap(const std::string name, std::ostream& ostr);

/// Read end notification to aid debugging.
void feasst_deserialize_endcap(const std::string name, std::istream& istr);

/// A factory method to construct objects from argtype
template <typename T>
std::shared_ptr<T> template_factory(
    std::map<std::string, std::shared_ptr<T> >& map,
    std::string class_name,
    argtype * args) {
  DEBUG("deserializing: " << class_name);
  if (map.count(class_name) == 0) {
    FATAL("The class name \"" << class_name << "\" is not recognized.");
  }
  std::shared_ptr<T> obj = map[class_name]->create(args);
  DEBUG("obj " << obj);
  return obj;
}

}  // namespace feasst

#endif  // FEASST_UTILS_SERIALIZE_H_
