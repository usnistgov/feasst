#include "utils/test/utils.h"
#include "monte_carlo/include/trial_factory.h"
#include "monte_carlo/include/trial_add.h"
#include "monte_carlo/include/trial_remove.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_rotate.h"

namespace feasst {

TEST(TrialFactory, serialize) {
  TrialFactory factory;
  factory.add(MakeTrialAdd());
  factory.add(MakeTrialRemove());
  factory.add(MakeTrialTranslate());
  factory.add(MakeTrialRotate());
//  std::stringstream ss;
//  factory.serialize(ss);
//  INFO(ss.str());
//  TrialFactory factory3(ss);
//  factory3.serialize(ss);
//  INFO(ss.str());
  TrialFactory factory2 = test_serialize(factory);
}

}  // namespace feasst
