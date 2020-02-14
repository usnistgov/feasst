
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
      const bool is_position_held = true) override;
  void before_select() override {
    Perturb::before_select();
    anywhere_.before_select();
  }
  void finalize(System * system) override;
  void revert(System * system) override;

  // serialize
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
