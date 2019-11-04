
#ifndef FEASST_MONTE_CARLO_PERTURB_REMOVE_H_
#define FEASST_MONTE_CARLO_PERTURB_REMOVE_H_

#include "monte_carlo/include/perturb_anywhere.h"

namespace feasst {

/**
  Remove a particle from the system.
 */
class PerturbRemove : public Perturb {
 public:
  PerturbRemove();

  void perturb(
      System * system,
      TrialSelect * select,
      Random * random,
      const bool is_position_held = true
      ) override {
    set_finalize_possible(true, select);

    if (is_position_held) {
      anywhere_.set_revert_possible(false, NULL);
    } else {
      anywhere_.perturb(system, select, random, is_position_held);
      set_revert_possible(true, select);
    }

    // setting trial state should go last so other perturbs do not overwrite
    select->set_trial_state(2);
  }

  void finalize(System * system) override {
    system->finalize();
    if (finalize_possible()) {
      system->get_configuration()->remove_particles(finalize_select()->mobile());
      // system->revert();
    }
  }

  void revert(System * system) override {
    if (revert_possible()) {
      anywhere_.revert(system);
    }
  }

  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbRemove(std::istream& istr);
  virtual ~PerturbRemove() {}

 private:
  // temporary
  PerturbAnywhere anywhere_;
};

inline std::shared_ptr<PerturbRemove> MakePerturbRemove() {
  return std::make_shared<PerturbRemove>();
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_REMOVE_H_
