#include "utils/test/utils.h"
#include "configuration/include/properties.h"
#include "utils/include/debug.h"
#include "math/include/constants.h"

namespace feasst {

TEST(Properties, properties) {
  Properties properties;
  TRY(
    properties.value("bananas");
    CATCH_PHRASE("not found");
  );
  double value;
  EXPECT_FALSE(properties.value("bananas", &value));
  properties.add_or_set("bananas", 12);
  EXPECT_TRUE(properties.value("bananas", &value));
  EXPECT_NEAR(properties.value("bananas"), 12, NEAR_ZERO);
  properties.add_or_set("bananas", 11);
  EXPECT_NEAR(properties.value("bananas"), 11, NEAR_ZERO);
  TRY(
    properties.set("apples", 2.3);
    CATCH_PHRASE("not found");
  );
  properties.set("bananas", 2.3);
  EXPECT_NEAR(properties.value("bananas"), 2.3, NEAR_ZERO);
  TRY(
    properties.add("bananas", 12);
    CATCH_PHRASE("already exists");
  );
  properties.check();

  // serialize
  properties.add("now", -10325231);
  Properties prop2 = test_serialize(properties);
  EXPECT_EQ(prop2.str(), properties.str());
}

}  // namespace feasst
