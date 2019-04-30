#include <gtest/gtest.h>
#include "example/include/model_example.h"

namespace feasst {

TEST(ModelExample, model_example) {
  ModelExample example;

  // test energy
  EXPECT_EQ(0., example.energy(0, 0, 0, ModelParams()));

  // test serialization
  std::stringstream ss, ss2;
  example.serialize(ss);
  std::shared_ptr<Model> model = example.deserialize(ss);
  model->serialize(ss2);
  EXPECT_EQ(ss.str(), ss2.str());
}

}  // namespace feasst
