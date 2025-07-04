
#ifndef FEASST_CHAIN_SELECT_TWO_SITES_H_
#define FEASST_CHAIN_SELECT_TWO_SITES_H_

#include <memory>
#include "monte_carlo/include/trial_select.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  A random particle of given type is selected if previously perturbed sites
  are not available.
  Select a single bond from given anchor to mobile sites.
 */
class SelectTwoSites : public TrialSelect {
 public:
  //@{
  /** @name Arguments
    - mobile_site: name of the mobile site.
    - mobile_site2: name of the second mobile site.
    - particle_type2: type of particle for second site. If -1, use the same
      particle (default: -1).
   */
  explicit SelectTwoSites(argtype args = argtype());
  explicit SelectTwoSites(argtype * arg);

  //@}
  /** @name Public Functions
   */
  //@{

  /// mobile is sized.
  void precompute(System * system) override;

  bool select(const Select& perturbed,
    System* system,
    Random * random,
    TrialSelect * previous_select) override;

  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit SelectTwoSites(std::istream& istr);
  virtual ~SelectTwoSites() {}

  //@}
 protected:
  void serialize_select_two_sites_(std::ostream& ostr) const;

 private:
  int particle_type2_;
  std::string mobile_site_name_, mobile_site2_name_;
};

inline std::shared_ptr<SelectTwoSites> MakeSelectTwoSites(
    argtype args = argtype()) {
  return std::make_shared<SelectTwoSites>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_SELECT_TWO_SITES_H_
