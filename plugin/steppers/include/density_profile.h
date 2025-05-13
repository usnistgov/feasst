
#ifndef FEASST_STEPPERS_DENSITY_PROFILE_H_
#define FEASST_STEPPERS_DENSITY_PROFILE_H_

#include "math/include/histogram.h"
#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  Compute the Density as a function of distance, r, in a given dimension for
  each site type.
 */
class DensityProfile : public Analyze {
 public:
  //@{
  /** @name Arguments
    - dimension: profile runs along this dimension (default: 0).
    - dr: profile bin size (default: 0.1).
    - center: specify the center of a single bin (default: 0.).
   */
  explicit DensityProfile(argtype args = argtype());
  explicit DensityProfile(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::string header(const MonteCarlo& mc) const override;
  void initialize(MonteCarlo * mc) override;
  void update(const MonteCarlo& mc) override;

  /// Return the profile for a given site type.
  /// The first index is the bin.
  /// The second index is the type.
  /// The third index is 0 - r, 1 - value.
  std::vector<std::vector<std::vector<double> > > profile() const;

  std::string write(const MonteCarlo& mc) override;

  // serialize
  std::string class_name() const override { return std::string("DensityProfile"); }
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<DensityProfile>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<DensityProfile>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit DensityProfile(std::istream& istr);
  explicit DensityProfile(const Analyze& density_profile);

  //@}
 private:
  int dimension_;
  double dr_, center_;
  std::vector<Histogram> data_;
};

inline std::shared_ptr<DensityProfile> MakeDensityProfile(
    argtype args = argtype()) {
  return std::make_shared<DensityProfile>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_DENSITY_PROFILE_H_
