
#ifndef FEASST_CONFINEMENT_HENRY_COEFFICIENT_H_
#define FEASST_CONFINEMENT_HENRY_COEFFICIENT_H_

#include "monte_carlo/include/analyze.h"
#include "math/include/accumulator.h"

namespace feasst {

/**
  Assumes there is only one Trial, TrialAdd, and it is AlwaysReject-ed.
  Accumulate \f$\langle e^{-\beta \Delta U}\rangle\f$,
  where \f$\Delta U\f$ is the energy contribution of the attempt to add the
  particle.
 */
class HenryCoefficient : public Analyze {
 public:
  //@{
  /** @name Arguments
    - num_beta_taylor: number of derivatives of second virial ratio with
      respect to beta. (default: 0).
   */
  explicit HenryCoefficient(argtype args = argtype());
  explicit HenryCoefficient(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Return the number of beta derivatives, starting with 1.
  int num_beta_taylor() const { return static_cast<int>(beta_taylor_.size()); }

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

  const Accumulator& coefficient() const { return accumulator(); }

  // serialize
  std::string class_name() const override { return std::string("HenryCoefficient"); }
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<HenryCoefficient>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit HenryCoefficient(std::istream& istr);
  virtual ~HenryCoefficient() {}
  //@}

 private:
  std::vector<Accumulator> beta_taylor_;
};

inline std::shared_ptr<HenryCoefficient> MakeHenryCoefficient(
    argtype args = argtype()) {
  return std::make_shared<HenryCoefficient>(args);
}

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_HENRY_COEFFICIENT_H_
