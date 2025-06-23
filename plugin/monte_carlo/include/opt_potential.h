
#ifndef FEASST_MONTE_CARLO_OPT_POTENTIAL_H_
#define FEASST_MONTE_CARLO_OPT_POTENTIAL_H_

#include <memory>
#include <vector>
#include "monte_carlo/include/action.h"

namespace feasst {

/**
  OptPotential does not optimize or speed up the potential.
  CheckEnergy simply ensures this optimized Potential returns the same energy as
  an unoptimized Potential.
  Thus, OptPotential can be used for testing alternative implementations which
  may be faster than the more simple and easier to test algorithms.
  For example, a cell list may be implemented in the OptPotential but not in the
  Potential to ensure that the cell list is not missing important interactions.
  See <a href="../../monte_carlo/tutorial/tutorial_4_lj_triclinic_celllist.html">this tutorial</a>
  for an example of testing VisitModelCell.
  If there are multiple Potential, OptPotential must include all Potential
  even if identical (e.g., in the example, LongRangeCorrections are required in
  OptPotential if included in Potential).
 */
class OptPotential : public Action {
 public:
  //@{
  /** @name Arguments
    - Potential arguments.
   */
  explicit OptPotential(argtype args = argtype());
  explicit OptPotential(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<OptPotential>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<OptPotential>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit OptPotential(std::istream& istr);
  virtual ~OptPotential() {}

  //@}
 private:
  int configuration_index_;
  argtype args_;
};

inline std::shared_ptr<OptPotential> MakeOptPotential(
    argtype args = argtype()) {
  return std::make_shared<OptPotential>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_OPT_POTENTIAL_H_
