

#ifndef FEASST_UTILS_SERIALIZE_EXTRA_H_
#define FEASST_UTILS_SERIALIZE_EXTRA_H_

#include <memory>
#include <string>
#include "utils/include/serialize.h"

namespace feasst {

/// Return a shared pointer to the base class of model after construction of
/// the full derived class.
/// see https://isocpp.org/wiki/faq/serialization
template <typename T>
std::shared_ptr<T> template_deserialize(
    std::map<std::string, std::shared_ptr<T> >& map,  // NOLINT
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

/// A factory method to construct objects from argtype
template <typename T>
std::shared_ptr<T> template_factory(
    std::map<std::string, std::shared_ptr<T> >& map,  // NOLINT
    std::string class_name,
    argtype * args) {
  DEBUG("deserializing: " << class_name);
  if (map.count(class_name) == 0) {
    INFO("candidates:");
    for (const auto& ele : map) {
      INFO(ele.first);
    }
    FATAL("The class name \"" << class_name << "\" is not recognized.");
  }
  std::shared_ptr<T> obj = map[class_name]->create(args);
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

#endif  // FEASST_UTILS_SERIALIZE_EXTRA_H_
