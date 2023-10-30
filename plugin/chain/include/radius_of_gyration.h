
#ifndef FEASST_CHAIN_RADIUS_OF_GYRATION
#define FEASST_CHAIN_RADIUS_OF_GYRATION

#include "math/include/accumulator.h"
#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  Accumulate average radius of gyration.
 */
class RadiusOfGyration : public Analyze {
 public:
  //@{
  /** @name Arguments
    - group_index: index of the Configuration::group (default: 0).
    - Stepper arguments.
   */
  explicit RadiusOfGyration(argtype args = argtype());
  explicit RadiusOfGyration(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::string header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trials) const override;

  void initialize(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override;

  void update(const Criteria& criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  std::string write(const Criteria& criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  /// Return the squared radius of gyration
  const Accumulator& radius_of_gyration() const { return accumulator(); }

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
};

inline std::shared_ptr<RadiusOfGyration> MakeRadiusOfGyration(argtype args = argtype()) {
  return std::make_shared<RadiusOfGyration>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_RADIUS_OF_GYRATION
