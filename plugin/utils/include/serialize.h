
#ifndef FEASST_UTILS_SERIALIZE_H_
#define FEASST_UTILS_SERIALIZE_H_

#include <string>
#include <sstream>
#include <vector>
#include <deque>
#include <memory>
#include <map>
#include "utils/include/debug.h"
#include "utils/include/utils_io.h"

/**
  Utility functions for serialization of objects into human-readable character
  streams.
  Each serialization function has an accompying deserialization function which
  performs the inverse operation.

  Note that some deserialization functions of derived objects do not work as
  the commented templates.
  But when the source code is copied to the specific object's deserialization
  function, then it does work, for unknown reasons.
 */

namespace feasst {

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

/// Serialize the 3D vector.
template <typename T>
void feasst_serialize(const std::vector<std::vector<std::vector<T> > >& vector,
    std::ostream& ostr) {
  ostr << MAX_PRECISION;
  ostr << vector.size() << " ";
  for (const std::vector<std::vector<T> >& vec2 : vector) {
    ostr << vec2.size() << " ";
    for (const std::vector<T>& vec3 : vec2) {
      ostr << vec3.size() << " ";
      for (const T& element : vec3) {
        ostr << element << " ";
      }
    }
  }
}

/// Deserialize the 3D vector.
template <typename T>
void feasst_deserialize(std::vector<std::vector<std::vector<T> > > * vector,
    std::istream& istr) {
  int dim1;
  istr >> dim1;
  vector->resize(dim1);
  for (int index1 = 0; index1 < dim1; ++index1) {
    int dim2;
    istr >> dim2;
    (*vector)[index1].resize(dim2);
    for (int index2 = 0; index2 < dim2; ++index2) {
      int dim3;
      istr >> dim3;
      (*vector)[index1][index2].resize(dim3);
      for (int index3 = 0; index3 < dim3; ++index3) {
        istr >> (*vector)[index1][index2][index3];
      }
    }
  }
}

/// Serialize the 4D vector.
template <typename T>
void feasst_serialize(const std::vector<std::vector<std::vector<std::vector<
    T> > > >& vector,
    std::ostream& ostr) {
  ostr << MAX_PRECISION;
  ostr << vector.size() << " ";
  for (const std::vector<std::vector<std::vector<T> > >& vec2 : vector) {
    ostr << vec2.size() << " ";
    for (const std::vector<std::vector<T> >& vec3 : vec2) {
      ostr << vec3.size() << " ";
      for (const std::vector<T>& vec4 : vec3) {
        ostr << vec4.size() << " ";
        for (const T& element : vec4) {
          ostr << element << " ";
        }
      }
    }
  }
}

/// Deserialize the 4D vector.
template <typename T>
void feasst_deserialize(std::vector<std::vector<std::vector<std::vector<
    T> > > > * vector,
    std::istream& istr) {
  int dim1;
  istr >> dim1;
  vector->resize(dim1);
  for (int index1 = 0; index1 < dim1; ++index1) {
    int dim2;
    istr >> dim2;
    (*vector)[index1].resize(dim2);
    for (int index2 = 0; index2 < dim2; ++index2) {
      int dim3;
      istr >> dim3;
      (*vector)[index1][index2].resize(dim3);
      for (int index3 = 0; index3 < dim3; ++index3) {
        int dim4;
        istr >> dim4;
        (*vector)[index1][index2][index3].resize(dim4);
        for (int index4 = 0; index4 < dim4; ++index4) {
          istr >> (*vector)[index1][index2][index3][index4];
        }
      }
    }
  }
}

/// Serialize the 5D vector.
template <typename T>
void feasst_serialize(const std::vector<std::vector<std::vector<std::vector<
    std::vector<T> > > > >& vector,
    std::ostream& ostr) {
  ostr << MAX_PRECISION;
  ostr << vector.size() << " ";
  for (const std::vector<std::vector<std::vector<std::vector<T> > > >& vec2
    : vector) {
    ostr << vec2.size() << " ";
    for (const std::vector<std::vector<std::vector<T> > >& vec3 : vec2) {
      ostr << vec3.size() << " ";
      for (const std::vector<std::vector<T> >& vec4 : vec3) {
        ostr << vec4.size() << " ";
        for (const std::vector<T>& vec5 : vec4) {
          ostr << vec5.size() << " ";
          for (const T& element : vec5) {
            ostr << element << " ";
          }
        }
      }
    }
  }
}

/// Deserialize the 5D vector.
template <typename T>
void feasst_deserialize(std::vector<std::vector<std::vector<std::vector<
    std::vector<T> > > > > * vector,
    std::istream& istr) {
  int dim1;
  istr >> dim1;
  vector->resize(dim1);
  for (int index1 = 0; index1 < dim1; ++index1) {
    int dim2;
    istr >> dim2;
    (*vector)[index1].resize(dim2);
    for (int index2 = 0; index2 < dim2; ++index2) {
      int dim3;
      istr >> dim3;
      (*vector)[index1][index2].resize(dim3);
      for (int index3 = 0; index3 < dim3; ++index3) {
        int dim4;
        istr >> dim4;
        (*vector)[index1][index2][index3].resize(dim4);
        for (int index4 = 0; index4 < dim4; ++index4) {
          int dim5;
          istr >> dim5;
          (*vector)[index1][index2][index3][index4].resize(dim5);
          for (int index5 = 0; index5 < dim5; ++index5) {
            istr >> (*vector)[index1][index2][index3][index4][index5];
          }
        }
      }
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
    << "Otherwise, this is likely due to the lack of a static mapper "
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

}  // namespace feasst

#endif  // FEASST_UTILS_SERIALIZE_H_
