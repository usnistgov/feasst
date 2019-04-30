#include <sstream>
#include <gtest/gtest.h>
#include "core/include/model_square_well.h"

namespace feasst {

TEST(ModelSquareWell, serialize) {
  ModelSquareWell model;
  std::stringstream istr;
  model.serialize(istr);
  EXPECT_EQ("ModelSquareWell 1 ", istr.str());
  std::shared_ptr<Model> model2 = ModelSquareWell().deserialize(istr);
  istr.str("");
  model2->serialize(istr);
  EXPECT_EQ("ModelSquareWell 1 ", istr.str());
}

}  // namespace feasst
