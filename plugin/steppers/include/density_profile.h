
#ifndef FEASST_STEPPERS_DENSITY_PROFILE_H_
#define FEASST_STEPPERS_DENSITY_PROFILE_H_

#include "math/include/histogram.h"
#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  Density profile
 */
class DensityProfile : public Analyze {
 public:
  /**
    args:
    - dimension: profile runs along this dimension (default: 0).
    - dr: profile bin size (default: 0.1).
    - center: specify the center of a single bin (default: 0.).
   */
  explicit DensityProfile(argtype args = argtype());

  std::string header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trials) const override;

  void initialize(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override;

  void update(const Criteria& criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  /// Return the profile for a given type.
  /// The first index is the bin.
  /// The second index is the type.
  /// The third index is 0 - r, 1 - value.
  std::vector<std::vector<std::vector<double> > > profile() const;

  std::string write(const Criteria& criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  // serialize
  std::string class_name() const override { return std::string("DensityProfile"); }
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<DensityProfile>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit DensityProfile(std::istream& istr);
  explicit DensityProfile(const Analyze& density_profile);

 private:
  int dimension_;
  double dr_, center_;

  // temporary and not serialized
  std::vector<Histogram> data_;
};

inline std::shared_ptr<DensityProfile> MakeDensityProfile(
    argtype args = argtype()) {
  return std::make_shared<DensityProfile>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_DENSITY_PROFILE_H_
