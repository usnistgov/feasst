
#ifndef FEASST_CHAIN_SELECT_PARTICLE_PIVOT_H_
#define FEASST_CHAIN_SELECT_PARTICLE_PIVOT_H_

#include <vector>
#include <memory>
#include "monte_carlo/include/trial_select.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  Random selection a whole particle of given type to pivot.
  Exclude the pivot site from the mobile selection for optimized energy calc.
 */
class SelectParticlePivot : public TrialSelect {
 public:
  //@{
  /** @name Arguments
    - pivot_site: the site name to use as the center of the rotation.
      By default, the pivot site is the first site listed in the particle.
   */
  explicit SelectParticlePivot(argtype args = argtype());
  explicit SelectParticlePivot(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void precompute(System * system) override;

  bool select(const Select& perturbed,
    System* system,
    Random * random,
    TrialSelect * previous_select) override;

  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit SelectParticlePivot(std::istream& istr);
  virtual ~SelectParticlePivot() {}

  //@}
 protected:
  void serialize_select_particle_pivot_(std::ostream& ostr) const;

 private:
  int pivot_site_ = 0;
  std::string pivot_site_name_;
};

inline std::shared_ptr<SelectParticlePivot> MakeSelectParticlePivot(
    argtype args = argtype()) {
  return std::make_shared<SelectParticlePivot>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_SELECT_PARTICLE_PIVOT_H_
