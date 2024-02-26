
#ifndef FEASST_CHAIN_PERTURB_REPTATE_H_
#define FEASST_CHAIN_PERTURB_REPTATE_H_

#include "monte_carlo/include/perturb_distance.h"

namespace feasst {

/**
  For a reptation, if new bond is accepted, then change the positions of all the
  sites along the chain.
  For heteropolymers, this perturbation changes the composition.
 */
class PerturbReptate : public PerturbDistance {
 public:
  PerturbReptate(argtype args = argtype()) : PerturbReptate(&args) {
    FEASST_CHECK_ALL_USED(args);
  }
  PerturbReptate(argtype * args) : PerturbDistance(args) {
    class_name_ = "PerturbReptate";
  }
  void move(const bool is_position_held, System * system, TrialSelect * select,
      Random * random, Acceptance * acceptance) override {
    PerturbDistance::move(is_position_held, system, select, random, acceptance);
    set_finalize_possible(true, select);
  }

  void finalize(System * system) override;
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbReptate(std::istream& istr);
  virtual ~PerturbReptate() {}

 protected:
  void serialize_perturb_reptate_(std::ostream& ostr) const;
};

inline std::shared_ptr<PerturbReptate> MakePerturbReptate(argtype args = argtype()) {
  return std::make_shared<PerturbReptate>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_PERTURB_REPTATE_H_
