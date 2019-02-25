#include <gtest/gtest.h>
#include "core/include/properties.h"
#include "core/include/debug.h"
#include "core/include/constants.h"

namespace feasst {

TEST(Properties, properties) {
  Properties properties;
  try {
    properties.value("bananas");
    CATCH_PHRASE("not found");
  }
  double value;
  EXPECT_FALSE(properties.value("bananas", &value));
  properties.add_or_set("bananas", 12);
  EXPECT_TRUE(properties.value("bananas", &value));
  EXPECT_NEAR(properties.value("bananas"), 12, NEAR_ZERO);
  properties.add_or_set("bananas", 11);
  EXPECT_NEAR(properties.value("bananas"), 11, NEAR_ZERO);
  try {
    properties.set("apples", 2.3);
    CATCH_PHRASE("not found");
  }
  properties.set("bananas", 2.3);
  EXPECT_NEAR(properties.value("bananas"), 2.3, NEAR_ZERO);
  try {
    properties.add("bananas", 12);
    CATCH_PHRASE("already exists");
  }
  properties.check_size();
}

}  // namespace feasst
