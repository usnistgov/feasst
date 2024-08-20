
#ifndef FEASST_STEPPERS_SCATTERING_H_
#define FEASST_STEPPERS_SCATTERING_H_

#include <vector>
#include "math/include/position.h"
#include "math/include/accumulator.h"
#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  Compute the scattering intensity.

  \f$F(\vec{q}) = \sum f_i(|\vec{q}|R_i) \exp(-i \vec{q} \cdot \vec{r}_i)\f$

  where \f$r\f$ is the position of site \f$i\f$, \f$q\f$ is the frequency and
  \f$f_i\f$ are the scattering lengths, assumed to be
  spherical with a radius, \f$R_i = \sigma_i\f$.

  \f$f(qR) = 3 V(R) \frac{\sin(qR)-qR\cos(qr)}{q^3R^3}\f$,

  \f$V(R)=\frac{4}{3}\pi R^3\f$.

  The intensity, \f$I\f$, is then given by the \f$F(\vec{q})\f$ multiplied by
  its complex conjugate:

  \f$I(q) = |F(\vec{q})|^2\f$.
 */
class Scattering : public Analyze {
 public:
  //@{
  /** @name Arguments
    - num_frequency: the number of linearly spaced frequencies between the
      largest and the smallest, 2*pi/minimum_domain_length (default: 100).
    - Stepper arguments.
  */
  explicit Scattering(argtype args = argtype());
  explicit Scattering(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

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
  std::string class_name() const override { return std::string("Scattering"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<Scattering>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<Scattering>(args); }
  explicit Scattering(std::istream& istr);

  //@}
 private:
  int num_frequency_;
  std::vector<Position> kvecs_;
  std::vector<std::vector<double> > site_ff_;
  std::vector<Accumulator> iq_;

//  std::vector<double> iq_() const;
};

inline std::shared_ptr<Scattering> MakeScattering(
    argtype args = argtype()) {
  return std::make_shared<Scattering>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_SCATTERING_H_
