
#ifndef FEASST_MONTE_CARLO_PERTURB_ADD_REMOVE_H_
#define FEASST_MONTE_CARLO_PERTURB_ADD_REMOVE_H_

#include <string>
#include <memory>
#include "monte_carlo/include/perturb.h"

namespace feasst {

class PerturbAdd;
class PerturbRemove;

/**
  Add or remove a particle with equal probability.
 */
class PerturbAddRemove : public Perturb {
 public:
  explicit PerturbAddRemove(argtype args = argtype());
  explicit PerturbAddRemove(argtype * args);

  void precompute(TrialSelect * select, System * system) override;
  void before_select() override;
  void perturb(System * system, TrialSelect * select, Random * random,
      const bool is_position_held = false,
      Acceptance * acceptance = NULL) override;
  void revert(System * system) override;
  void finalize(System * system) override;
  std::string status_header() const override;
  std::string status() const override;

  // serialize
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbAddRemove(std::istream& istr);
  virtual ~PerturbAddRemove();

 private:
  std::unique_ptr<PerturbAdd> add_;
  std::unique_ptr<PerturbRemove> rm_;
  bool adding_;
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_ADD_REMOVE_H_
