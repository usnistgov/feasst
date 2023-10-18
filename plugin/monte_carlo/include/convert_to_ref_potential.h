
#ifndef FEASST_MONTE_CARLO_CONVERT_TO_REF_POTENTIAL_H_
#define FEASST_MONTE_CARLO_CONVERT_TO_REF_POTENTIAL_H_

#include <vector>
#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/action.h"

namespace feasst {

/**
  Make a new reference potential based on an existing potential.
 */
class ConvertToRefPotential : public Action {
 public:
  //@{
  /** @name Arguments
   */

  /**
    args:
    - potential_index: index of the full potential to copy as a template
      (default: 0).
    - cutoff: set cutoff of all site_types to this value.
      Ignore if -1 (default: -1).
    - use_cell: use VisitModelCell with cutoff as min_length (default: false).
   */
  explicit ConvertToRefPotential(argtype args = argtype());
  explicit ConvertToRefPotential(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<ConvertToRefPotential>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<ConvertToRefPotential>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit ConvertToRefPotential(std::istream& istr);
  virtual ~ConvertToRefPotential() {}

  //@}
 private:
  int potential_index_;
  double cutoff_;
  bool use_cell_;
};

inline std::shared_ptr<ConvertToRefPotential> MakeConvertToRefPotential(
    argtype args = argtype()) {
  return std::make_shared<ConvertToRefPotential>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_CONVERT_TO_REF_POTENTIAL_H_
