
#ifndef FEASST_MONTE_CARLO_PERTURB_ANYWHERE_H_
#define FEASST_MONTE_CARLO_PERTURB_ANYWHERE_H_

#include <memory>
#include "monte_carlo/include/perturb_translate.h"
#include "monte_carlo/include/perturb_rotate.h"

namespace feasst {

/// Rigidly move rigidly anywhere in the box with any orientation.
class PerturbAnywhere : public PerturbMove {
 public:
  PerturbAnywhere();

  // assumes first particle center is placed on origin before translation to center
  void set_position(const Position& center, System * system, TrialSelect * select) {
    Position add = center;
    add.subtract(select->mobile().particle_positions()[0]);
    translate_.move(add, system, select);
  }

  void move(System * system, TrialSelect * select, Random * random) override {
    ASSERT(std::abs(rotate_.tunable().value() - 180.) < NEAR_ZERO,
      "rotation tunable should be 180");
    rotate_.move(system, select, random);
    system->configuration().domain().random_position(&random_in_box_, random);
    set_position(random_in_box_, system, select);
    DEBUG("anywhere: " << random_in_box_.str());
  }

  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbAnywhere(std::istream& istr);
  virtual ~PerturbAnywhere() {}

 private:
  PerturbRotate rotate_;

  // temporary
  PerturbTranslate translate_;
  Position random_in_box_;
};

inline std::shared_ptr<PerturbAnywhere> MakePerturbAnywhere() {
  return std::make_shared<PerturbAnywhere>();
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_ANYWHERE_H_
