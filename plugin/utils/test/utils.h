#include <cmath>
#include <string>
#include <memory>
#include <gtest/gtest.h>
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "utils/include/io.h"
#include "utils/include/arguments.h"
#include "math/include/position.h"
#include "math/include/constants.h"

namespace feasst {

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
  DEBUG(ss.str());
  if (!expected.empty()) {
    EXPECT_EQ(ss.str(), expected);
  }
  auto object2 = T(ss);
  object2.serialize(ss2);
  DEBUG(ss.str());
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
  //INFO("ss " << ss.str());
  //INFO("ss2 " << ss2.str());
  if (compare) {
    EXPECT_EQ(ss.str(), ss2.str());
  }
  return object2;
}

/// Return a copy of the object via serialization.
/// Test that the serialization of the new object is the same.
template <typename T>
std::unique_ptr<T> test_serialize(const std::unique_ptr<T>& object,
    /// If not empty, check that the serialization matches expectation.
    std::string expected = "",
    /// Set to true to compare serialization of new object with old.
    bool compare = true) {
  std::stringstream ss, ss2;
  object->serialize(ss);
  DEBUG(ss.str());
  if (!expected.empty()) {
    EXPECT_EQ(ss.str(), expected);
  }
  std::unique_ptr<T> object2 = std::make_unique<T>(ss);
  object2->serialize(ss2);
  DEBUG(ss.str());
  if (compare) {
    EXPECT_EQ(ss.str(), ss2.str());
  }
  return object2;
}

/// Return a copy of the object via serialization.
/// Test that the serialization of the new object is the same.
template <typename T>
std::unique_ptr<T> test_serialize_unique(const T& object,
    /// If not empty, check that the serialization matches expectation.
    std::string expected = "",
    /// Set to true to compare serialization of new object with old.
    bool compare = true) {
  std::stringstream ss, ss2;
  object.serialize(ss);
  DEBUG(ss.str());
  if (!expected.empty()) {
    EXPECT_EQ(ss.str(), expected);
  }
  std::unique_ptr<T> object2 = std::make_unique<T>(ss);
  object2->serialize(ss2);
  DEBUG(ss.str());
  if (compare) {
    EXPECT_EQ(ss.str(), ss2.str());
  }
  return object2;
}

}  // namespace feasst
