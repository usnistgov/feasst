
#ifndef FEASST_STEPPERS_SCATTERING_H_
#define FEASST_STEPPERS_SCATTERING_H_

#include <vector>
#include "math/include/position.h"
#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  Compute the scattering intensity using FFTW (fftw.org) by gridding the system.
  Enabled for anisotropic interactions using VisitModelInnerTable.
 */
class ScatteringFFTW : public Analyze {
 public:
  /**
    args:
    - num_frequency: the number of linearly spaced frequencies between the
      largest and the smallest, 2*pi/minimum_domain_length (default: 100).
  */
  explicit ScatteringFFTW(argtype args = argtype());
  explicit ScatteringFFTW(argtype * args);

  void initialize(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override;

  void update(const Criteria& criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  std::string write(const Criteria& criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  int num_vectors() const { return static_cast<int>(kvecs_.size()); }

  // serialize
  std::string class_name() const override { return std::string("ScatteringFFTW"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<ScatteringFFTW>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<ScatteringFFTW>(args); }
  explicit ScatteringFFTW(std::istream& istr);

 private:
  int num_frequency_;
  std::vector<Position> kvecs_;
  std::vector<std::vector<double> > site_ff_;
  std::vector<Accumulator> iq_;

//  std::vector<double> iq_() const;
};

inline std::shared_ptr<ScatteringFFTW> MakeScatteringFFTW(
    argtype args = argtype()) {
  return std::make_shared<ScatteringFFTW>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_SCATTERING_H_
