#include <string>
#include <memory>
#include <gtest/gtest.h>
#include "utils/include/debug.h"

/// Return a copy of the object via serialization.
/// Test that the serialization of the new object is the same.
template <typename T>
T test_serialize(const T& object,
    /// If not empty, check that the serialization matches expectation.
    std::string expected = "",
    /// Set to true to compare serialization of new object with old.
    bool compare = true) {
  std::stringstream ss, ss2;
  object.serialize(ss);
  if (!expected.empty()) {
    EXPECT_EQ(ss.str(), expected);
  }
  auto object2 = T(ss);
  object2.serialize(ss2);
  // INFO(ss.str());
  if (compare) {
    EXPECT_EQ(ss.str(), ss2.str());
  }
  return object2;
}

template <typename T>
T test_serialize_no_comp(const T& object) {
  return test_serialize(object, "", false);
}

/// Return a copy of the object via serialization.
/// Test that the serialization of the new object is the same.
template <class T1, class T2>
std::shared_ptr<T2> test_serialize(T1 object,
    /// If not empty, check that the serialization matches expectation.
    std::string expected = "",
    /// Set to true to compare serialization of new object with old.
    bool compare = true) {
  std::stringstream ss, ss2;
  object.serialize(ss);
  if (!expected.empty()) {
    EXPECT_EQ(ss.str(), expected);
  }
  std::shared_ptr<T2> object2 = object.deserialize(ss);
  object2->serialize(ss2);
  // INFO(ss.str());
  if (compare) {
    EXPECT_EQ(ss.str(), ss2.str());
  }
  return object2;
}
