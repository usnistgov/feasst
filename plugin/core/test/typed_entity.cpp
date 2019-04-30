#include <gtest/gtest.h>
#include "core/include/typed_entity.h"

namespace feasst {

TEST(TypedEntity, serialize) {
  TypedEntity en;
  en.set_type(123);

  // serialize
  std::stringstream ss;
  en.serialize(ss);
  TypedEntity en2(ss);
  EXPECT_EQ(en.type(), 123);
  EXPECT_EQ(en.type(), en2.type());
}

}  // namespace feasst
