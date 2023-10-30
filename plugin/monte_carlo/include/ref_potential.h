
#ifndef FEASST_MONTE_CARLO_REF_POTENTIAL_H_
#define FEASST_MONTE_CARLO_REF_POTENTIAL_H_

#include <vector>
#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/action.h"

namespace feasst {

/**
  Add a reference potential.
  Unlike Potential, multiple RefPotential can use the same reference_index,
  which is 0 by default (see args).
 */
class RefPotential : public Action {
 public:
  //@{
  /** @name Arguments
    - reference_index: index of reference potential (default: 0).
    - configuration_index: index of configuration potential (default: 0).
    - Potential arguments.
   */
  explicit RefPotential(argtype args = argtype());
  explicit RefPotential(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<RefPotential>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<RefPotential>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit RefPotential(std::istream& istr);
  virtual ~RefPotential() {}

  //@}
 private:
  int reference_index_;
  int configuration_index_;
  argtype args_;
};

inline std::shared_ptr<RefPotential> MakeRefPotential(
    argtype args = argtype()) {
  return std::make_shared<RefPotential>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_REF_POTENTIAL_H_
