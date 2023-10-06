
#ifndef FEASST_GIBBS_PRESSURE_FROM_TEST_VOLUME_H_
#define FEASST_GIBBS_PRESSURE_FROM_TEST_VOLUME_H_

#include <vector>
#include <memory>
#include "utils/include/arguments.h"
#include "math/include/accumulator.h"
#include "math/include/histogram.h"
#include "monte_carlo/include/modify.h"

namespace feasst {

/**
  Compute pressure using test volume changes.
  See https://doi.org/10.1063/1.472721

  Hard systems must consider only volume reductions
  (i.e., negative delta_volume).

  \f$p = \frac{1}{\beta \Delta V}\ln \left\langle \left(\frac{V'}{V}\right)^N e^{-\beta\Delta U}\right\rangle \f$

  The output block standard deviation used the error propagation formula:

  \f$ f = \ln(g) \f$

  \f$ \sigma_f = \sigma_g/|g| \f$

  where \f$\sigma_g\f$ is a block average from the term above in the ensemble average.
 */
class PressureFromTestVolume : public Modify {
 public:
  /**
    args:
    - delta_volume: test volume change (default: 1e-4).
   */
  explicit PressureFromTestVolume(argtype args = argtype());
  explicit PressureFromTestVolume(argtype * args);

  std::string header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trials) const override;

  void initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) override;

  void update(Criteria * criteria,
    System * system,
    Random * random,
    TrialFactory * trial_factory) override;

  std::string write(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) override;

  // serialize
  std::string class_name() const override { return std::string("PressureFromTestVolume"); }
  std::shared_ptr<Modify> create(std::istream& istr) const override {
    return std::make_shared<PressureFromTestVolume>(istr); }
  std::shared_ptr<Modify> create(argtype * args) const override {
    return std::make_shared<PressureFromTestVolume>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit PressureFromTestVolume(std::istream& istr);

 private:
  double delta_volume_;
  Accumulator term_;
};

inline std::shared_ptr<PressureFromTestVolume> MakePressureFromTestVolume(
    argtype args = argtype()) {
  return std::make_shared<PressureFromTestVolume>(args);
}

}  // namespace feasst

#endif  // FEASST_GIBBS_PRESSURE_FROM_TEST_VOLUME_H_
