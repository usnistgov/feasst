
#ifndef FEASST_CHAIN_END_TO_END_DISTANCE
#define FEASST_CHAIN_END_TO_END_DISTANCE

#include "math/include/accumulator.h"
#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  Accumulate average end-to-end distance.
 */
class EndToEndDistance : public Analyze {
 public:
  //@{
  /** @name Arguments
   */

  /**
    args:
    - group_index: index of the Configuration::group (default: 0).
    - Stepper arguments.
   */
  explicit EndToEndDistance(argtype args = argtype());
  explicit EndToEndDistance(argtype * args);

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

  /// Return the end to end distance
  const Accumulator& end_to_end_distance() const { return accumulator(); }

  // serialize
  std::string class_name() const override { return std::string("EndToEndDistance"); }
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<EndToEndDistance>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<EndToEndDistance>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit EndToEndDistance(std::istream& istr);
  explicit EndToEndDistance(const Analyze& energy);

  //@}
 private:
  int group_index_;
};

inline std::shared_ptr<EndToEndDistance> MakeEndToEndDistance(argtype args = argtype()) {
  return std::make_shared<EndToEndDistance>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_END_TO_END_DISTANCE
