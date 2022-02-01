#include "utils/test/utils.h"
#include "monte_carlo/include/trial_translate.h"

namespace feasst {

TEST(TrialTranslate, serialize) {
  auto trial = MakeTrialTranslate();
//  std::stringstream ss;
//  trial.serialize(ss);
//  INFO(ss.str());
//  TrialTranslate trial2(ss);
//  trial2.serialize(ss);
//  INFO(ss.str());
  std::shared_ptr<Trial> trial2 = test_serialize<TrialTranslate, Trial>(*trial);
}

}  // namespace feasst
