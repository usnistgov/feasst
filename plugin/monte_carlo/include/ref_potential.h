
#ifndef FEASST_MONTE_CARLO_REF_POTENTIAL_H_
#define FEASST_MONTE_CARLO_REF_POTENTIAL_H_

#include <memory>
#include <vector>
#include "monte_carlo/include/action.h"

namespace feasst {

/**
  Add a reference potential.
  Each unique name added the potential to a different RefPotential that can be
  called specifically.
 */
class RefPotential : public Action {
 public:
  //@{
  /** @name Arguments
    - ref: name of reference potential (default: 0).
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
  std::string ref_;
  int configuration_index_;
  std::string config_;
  argtype args_;
};

inline std::shared_ptr<RefPotential> MakeRefPotential(
    argtype args = argtype()) {
  return std::make_shared<RefPotential>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_REF_POTENTIAL_H_
