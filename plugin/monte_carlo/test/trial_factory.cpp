#include "utils/test/utils.h"
#include "monte_carlo/include/trial_factory.h"
#include "monte_carlo/include/trial_rotate.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_add.h"
#include "monte_carlo/include/trial_remove.h"

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
  auto factory2 = test_serialize_unique(factory);
  TRY(
    factory2->remove(100);
    CATCH_PHRASE("Trial index:100 >= number of trials:4");
  );
}

}  // namespace feasst
