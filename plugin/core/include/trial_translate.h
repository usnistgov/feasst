
#ifndef FEASST_CORE_TRIAL_TRANSLATE_H_
#define FEASST_CORE_TRIAL_TRANSLATE_H_

#include "core/include/trial_move.h"
#include "core/include/perturb_translate.h"
#include "core/include/random.h"
#include "core/include/utils_io.h"

namespace feasst {

class TrialTranslate : public TrialMove {
 public:
  TrialTranslate() {
    perturb_trans_ = std::make_shared<PerturbTranslate>();
    set_perturb(perturb_trans_);
  }

  void move(System * system) override {
    const Position trajectory = random_.position_in_cube(
      system->dimension(),
      max_move()
    );
    perturb_trans_->translate_selection(trajectory, system);
  }

  virtual ~TrialTranslate() {}

 private:
  std::shared_ptr<PerturbTranslate> perturb_trans_;
  Random random_;
};

inline std::shared_ptr<TrialTranslate> TrialTranslateShrPtr() {
  return std::make_shared<TrialTranslate>();
}

}  // namespace feasst

#endif  // FEASST_CORE_TRIAL_TRANSLATE_H_
