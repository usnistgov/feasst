
#ifndef FEASST_UTILS_UTILS_IO_H_
#define FEASST_UTILS_UTILS_IO_H_

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <deque>
#include <numeric>
#include <memory>
#include <map>
#include "utils/include/debug.h"
#include "utils/include/utils.h"

using std::cout;
using std::endl;

namespace feasst {

/// Return string representation of vector
template<class T>
std::string feasst_str(const std::vector<T> &vec,
  /// use maximum precision
  const bool max_precision = false) {
  std::stringstream ss;
  if (max_precision) ss << MAX_PRECISION;
  for (unsigned int i = 0; i < vec.size(); ++i) {
    ss << vec[i] << " ";
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

/// Serialize bool.
void feasst_serialize(const bool val, std::ostream& ostr);

/// Deserialize bool.
void feasst_deserialize(bool * val, std::istream& istr);

/// Serialize string.
void feasst_serialize(const std::string str, std::ostream& ostr);

/// Deserialize string.
void feasst_deserialize(std::string * str, std::istream& istr);

/// Serialize generic to full precision.
template <typename T>
void feasst_serialize(const T val, std::ostream& ostr) {
  ostr << std::setprecision(std::numeric_limits<T>::digits10+2)
       << val << " ";
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

/// Serialize the 1D vector.
template <typename T>
void feasst_serialize(const std::vector<T>& vector, std::ostream& ostr) {
  ostr << MAX_PRECISION;
  ostr << vector.size() << " ";
  for (const T& element : vector) {
    ostr << element << " ";
  }
}

/// Deserialize the 1D vector.
template <typename T>
void feasst_deserialize(std::vector<T> * vector, std::istream& istr) {
  int num;
  istr >> num;
  vector->resize(num);
  for (int index = 0; index < num; ++index) {
    istr >> (*vector)[index];
  }
}

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

/// Serialize the 2D vector.
template <typename T>
void feasst_serialize(const std::vector<std::vector<T> >& vector,
    std::ostream& ostr) {
  ostr << MAX_PRECISION;
  ostr << vector.size() << " ";
  for (const std::vector<T>& inner_vector : vector) {
    ostr << inner_vector.size() << " ";
    for (const T& element : inner_vector) {
      ostr << element << " ";
    }
  }
}

/// Deserialize the 2D vector.
template <typename T>
void feasst_deserialize(std::vector<std::vector<T> > * vector,
    std::istream& istr) {
  int dim1;
  istr >> dim1;
  vector->resize(dim1);
  for (int index1 = 0; index1 < dim1; ++index1) {
    int dim2;
    istr >> dim2;
    (*vector)[index1].resize(dim2);
    for (int index2 = 0; index2 < dim2; ++index2) {
      istr >> (*vector)[index1][index2];
    }
  }
}

/// Serialize feasst object
template <typename T>
void feasst_serialize_fstobj(const T obj, std::ostream& ostr) {
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
///// Deserialize feasst derived object stored as shared pointer
//template <typename T>
//void feasst_deserialize_fstdr(std::shared_ptr<T> ptr, std::istream& istr) {
//  int existing;
//  istr >> existing;
//  if (existing != 0) {
//    ptr = ptr->deserialize(istr);
//  }
//}

/// Serialize vector of shared pointers of feasst objects
template <typename T>
void feasst_serialize(const std::vector<std::shared_ptr<T> >& vector,
    std::ostream& ostr) {
  ostr << vector.size() << " ";
  for (const std::shared_ptr<T> element : vector) {
    feasst_serialize(element, ostr);
  }
}

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

// HWH for unknown reasons, this function template does not work.
///// Deserialize vector of shared pointers of feasst derived objects
//template <typename T>
//void feasst_deserialize_fstdr(std::vector<std::shared_ptr<T> > * vector,
//    std::istream& istr) {
//  int dim1;
//  istr >> dim1;
//  vector->resize(dim1);
//  for (int index = 0; index < dim1; ++index) {
//    feasst_deserialize_fstdr((*vector)[index], istr);
//  }
//}

/// Return a shared pointer to the base class of model after construction of
/// the full derived class.
/// see https://isocpp.org/wiki/faq/serialization
template <typename T>
std::shared_ptr<T> template_deserialize(std::map<std::string, std::shared_ptr<T> >& map,
    std::istream& istr,
    /// Rewind istr position to read class name again (default: false).
    bool rewind = false) {
  std::string class_name;
  int pos = istr.tellg(); // record position before reading
  istr >> class_name;     // read class name

  // rewind position so constructors can reread class name.
  if (rewind) {
    istr.seekg(pos, istr.beg); // rewind to before reading the class name.
  }
  DEBUG("deserializing: " << class_name << " rewind? " << rewind);
  ASSERT(map.count(class_name) != 0, "The class name \"" << class_name << "\" "
    << "is not recognized during deserialization. "
    << "If the above class name is empty, there was a mis-match in stream. "
    << "Otherwise, this is likely due to the lack of a static mapper "
    << "which is typically implemented within the cpp file. "
    << "In rare cases, the absence of a constructor implementation inside "
    << "the cpp file possibly leads optimization to ignore the mapper.");
  std::shared_ptr<T> obj = map[class_name]->create(istr);
  DEBUG("obj " << obj);
  return obj;
  //return map[class_name]->create(istr);
}

/// Return a deep copy of a feasst derived class object.
/// This is implemented via serialization/deserialization.
template <typename T>
std::shared_ptr<T> deep_copy_derived(std::shared_ptr<T> object) {
  std::stringstream ss;
  object->serialize(ss);
  return object->deserialize(ss);
}

}  // namespace feasst

#endif  // FEASST_UTILS_UTILS_IO_H_
