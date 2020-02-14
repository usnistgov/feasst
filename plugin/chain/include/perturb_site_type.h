
#ifndef FEASST_MONTE_CARLO_PERTURB_SITE_TYPE_H_
#define FEASST_MONTE_CARLO_PERTURB_SITE_TYPE_H_

#include "monte_carlo/include/perturb_anywhere.h"

namespace feasst {

/**
  Change the type of a site.
 */
class PerturbSiteType : public Perturb {
 public:
  /**
    args:
    - type: type to set for site.
   */
  explicit PerturbSiteType(const argtype& args = argtype());

  //initialize ghost selection in TrialSelect?
  void precompute(TrialSelect * select, System * system) override {
    select->set_ghost(true); }

  /// Set the site type.
  void set_site_type(
    System * system,
    const Select& select,
    const int type);

  void perturb(
      System * system,
      TrialSelect * select,
      Random * random,
      const bool is_position_held = false) override;

  void revert(System * system) override;
  void finalize(System * system) override;
  std::string status_header() const override;
  std::string status() const override;

  // serialize
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbSiteType(std::istream& istr);
  virtual ~PerturbSiteType() {}

 private:
  int new_site_type_;

  // temporary
  int old_site_type_;
  int old_particle_type_;
};

inline std::shared_ptr<PerturbSiteType> MakePerturbSiteType(const argtype& args = argtype()) {
  return std::make_shared<PerturbSiteType>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_SITE_TYPE_H_
