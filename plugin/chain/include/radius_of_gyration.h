
#ifndef FEASST_CHAIN_RADIUS_OF_GYRATION
#define FEASST_CHAIN_RADIUS_OF_GYRATION

#include "math/include/accumulator.h"
#include "math/include/histogram.h"
#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  Accumulate average radius of gyration assuming unit mass for all sites.
 */
class RadiusOfGyration : public Analyze {
 public:
  //@{
  /** @name Arguments
    - group_index: index of the Configuration::group (default: 0).
    - print_histogram: if true, print histogram (default: false).
    - Histogram arguments.
    - Stepper arguments.
   */
  explicit RadiusOfGyration(argtype args = argtype());
  explicit RadiusOfGyration(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::string header(const MonteCarlo& mc) const override;
  void initialize(MonteCarlo * mc) override;
  void update(const MonteCarlo& mc) override;
  std::string write(const MonteCarlo& mc) override;

  /// Return the radius of gyration
  const Accumulator& radius_of_gyration() const { return accumulator(); }

  /// Return the histogram of the radius of gyration
  const std::shared_ptr<Histogram> histogram() const { return hist_; }

  /// Return the accumulator for radius of gyration times the energy for extrapolation
  const Accumulator& rg_e() const { return rg_e_; }
  const Accumulator& rg_e2() const { return rg_e2_; }

  // serialize
  std::string class_name() const override { return std::string("RadiusOfGyration"); }
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<RadiusOfGyration>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<RadiusOfGyration>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit RadiusOfGyration(std::istream& istr);
  explicit RadiusOfGyration(const Analyze& energy);

  //@}
 private:
  int group_index_;
  Accumulator rg_e_, rg_e2_;
  std::shared_ptr<Histogram> hist_;
};

inline std::shared_ptr<RadiusOfGyration> MakeRadiusOfGyration(argtype args = argtype()) {
  return std::make_shared<RadiusOfGyration>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_RADIUS_OF_GYRATION
