
#ifndef FEASST_MONTE_CARLO_CONSTRAIN_VOLUME_BY_CUTOFF_H_
#define FEASST_MONTE_CARLO_CONSTRAIN_VOLUME_BY_CUTOFF_H_

#include <memory>
#include "monte_carlo/include/constraint.h"

namespace feasst {

/**
  Constrain the Domain volume to be atleast twice the cutoff.
 */
class ConstrainVolumeByCutoff : public Constraint {
 public:
  //@{
  /** @name Arguments
   */
  explicit ConstrainVolumeByCutoff(argtype args = argtype());
  explicit ConstrainVolumeByCutoff(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  bool is_allowed(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) const override;

  std::shared_ptr<Constraint> create(std::istream& istr) const override {
    return std::make_shared<ConstrainVolumeByCutoff>(istr); }
  std::shared_ptr<Constraint> create(argtype * args) const override {
    return std::make_shared<ConstrainVolumeByCutoff>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit ConstrainVolumeByCutoff(std::istream& istr);
  virtual ~ConstrainVolumeByCutoff() {}

  //@}
 protected:
  void serialize_constrain_volume_by_cutoff_(std::ostream& ostr) const;
};

inline std::shared_ptr<ConstrainVolumeByCutoff> MakeConstrainVolumeByCutoff(
    argtype args = argtype()) {
  return std::make_shared<ConstrainVolumeByCutoff>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_CONSTRAIN_VOLUME_BY_CUTOFF_H_
