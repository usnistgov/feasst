
#ifndef FEASST_CONFINEMENT_HENRY_COEFFICIENT_H_
#define FEASST_CONFINEMENT_HENRY_COEFFICIENT_H_

#include "monte_carlo/include/analyze.h"
#include "math/include/accumulator.h"

namespace feasst {

/**
  Attempts TrialAdd, where all attempts are AlwaysReject ed.
  Assumes there is only one Trial, TrialAdd.
  Accumulates \f$\langle e^{-\beta \Delta U}\rangle\f$,
  where \f$\Delta U\f$ is the energy contribution of the attempt to add the
  particle.

  The purpose of this Analyze class is to enable computation of the Henry's Law constant for an adsorbate
  particle in a static medium. (The medium could be a porous material, hard confinment,
  position-dependent potential, frozen configuration of other particles, etc.)

  Optionally, when the num_beta_taylor argument is specified (see below), HenryCoefficient
  also measures the averages:

  \f$ K_j = \frac{\langle (-\Delta U)^j ~ e^{-\beta \Delta U}\rangle}{j!} \f$

  \f$ Q_j = \frac{\langle (-\Delta U)^j ~ e^{-2 \beta \Delta U}\rangle}{j! ~ j!} \f$

\rst
These quantities may be used to extrapolate the HenryCoefficient in beta space using a Taylor-series expansion.
See :footcite:t:`siderius_temperature_2022` for background information and derivation of the listed terms.

When num_beta_taylor is specified, the output file of HenryCoefficient includes a JSON-formatted
dictionary in the file header that stores relevant output:

- `beta_taylor` contains :math:`K_j` as a list for :math:`0 \leq j \leq` num_beta_taylor.
- `beta_taylor2` contains :math:`Q_j` as a list for :math:`0 \leq j \leq` 2*num_beta_taylor.
- `beta` stores the inverse temperature.
- `num_trials` records the number of trial insertions.

:math:`Q_j` is used for computation of the covariance matrix of the extrapolation coefficients, e.g.,

:math:`\textrm{cov} \left( K_j, K_l \right) = \frac{N_{mc}}{N_{mc}-1} \left( Q_{j+l} \frac{(j+l)! ~ (j+l)!}{j! ~ l!} - K_j \cdot K_l  \right)`

where :math:`N_{mc}` is the number of test insertions (equal to `num_trials` in the output dictionary).
The covariance matrix enable estimation of the uncertainty in an extrapolation function built from the Taylor-series expansion.

References:

.. footbibliography::
\endrst
 */
class HenryCoefficient : public Analyze {
 public:
  //@{
  /** @name Arguments
    - num_beta_taylor: number of derivatives of second virial ratio with
      respect to beta. (default: 0).
    - write_precision: number of decimals in writing taylor coefficients
      (default: 8).
    - Stepper arguments.
   */
  explicit HenryCoefficient(argtype args = argtype());
  explicit HenryCoefficient(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Return the number of beta derivatives, starting with 1.
  int num_beta_taylor() const { return static_cast<int>(beta_taylor_.size()); }

  std::string header(const MonteCarlo& mc) const override;
  void initialize(MonteCarlo * mc) override;
  void update(const MonteCarlo& mc) override;
  std::string write(const MonteCarlo& mc) override;

  const Accumulator& coefficient() const { return accumulator(); }

  // serialize
  std::string class_name() const override { return std::string("HenryCoefficient"); }
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<HenryCoefficient>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<HenryCoefficient>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit HenryCoefficient(std::istream& istr);
  virtual ~HenryCoefficient() {}
  //@}

 private:
  std::vector<Accumulator> beta_taylor_;
  std::vector<Accumulator> beta_taylor2_;
  int write_precision_;
};

inline std::shared_ptr<HenryCoefficient> MakeHenryCoefficient(
    argtype args = argtype()) {
  return std::make_shared<HenryCoefficient>(args);
}

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_HENRY_COEFFICIENT_H_
