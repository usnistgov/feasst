#include "utils/test/utils.h"
#include "configuration/include/typed_entity.h"

namespace feasst {

TEST(TypedEntity, serialize) {
  TypedEntity en;
  en.set_type(123);

  // serialize
  TypedEntity en2 = test_serialize(en);
  EXPECT_EQ(en.type(), 123);
  EXPECT_EQ(en.type(), en2.type());
}

}  // namespace feasst
