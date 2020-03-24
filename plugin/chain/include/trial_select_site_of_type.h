
#ifndef FEASST_CHAIN_TRIAL_SELECT_SITE_OF_TYPE_H_
#define FEASST_CHAIN_TRIAL_SELECT_SITE_OF_TYPE_H_

#include "monte_carlo/include/trial_select.h"

namespace feasst {

/// Select a random site of a given type.
class TrialSelectSiteOfType : public TrialSelect {
 public:
  /**
    args:
    - site_type: type of site to select.
   */
  TrialSelectSiteOfType(const argtype& args);

  int site_type() const { return site_type_; }

  /// Select a random site of given type in randomly selected particle.
  int random_site_in_particle(
      const Configuration& config,
      Select * select,
      Random * random);

  bool select(const Select& perturbed,
    System* system,
    Random * random) override;

  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialSelectSiteOfType(std::istream& istr);
  virtual ~TrialSelectSiteOfType() {}

 protected:
  void serialize_trial_select_segment_(std::ostream& ostr) const;

 private:
  int site_type_;
};

inline std::shared_ptr<TrialSelectSiteOfType> MakeTrialSelectSiteOfType(
    const argtype &args = argtype()) {
  return std::make_shared<TrialSelectSiteOfType>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_SELECT_SITE_OF_TYPE_H_
