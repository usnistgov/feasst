#include <sstream>
#include <gtest/gtest.h>
#include "core/include/model_hard_sphere.h"

namespace feasst {

TEST(ModelHardSphere, serialize) {
  ModelHardSphere model;
  std::stringstream ss, ss2;
  model.serialize(ss);
  EXPECT_EQ("ModelHardSphere 607 ", ss.str());
  std::shared_ptr<Model> model2 = ModelHardSphere().deserialize(ss);
  model2->serialize(ss2);
  EXPECT_EQ(ss.str(), ss2.str());
}

}  // namespace feasst
