
#ifndef FEASST_SHAPE_FORMULA_SINE_WAVE_H_
#define FEASST_SHAPE_FORMULA_SINE_WAVE_H_

#include <string>
#include <memory>
#include <vector>
#include "math/include/formula.h"

namespace feasst {

/**
  \f$y = shift + amplitude*\sin(2\pi (x - phase)/ width)\f$
 */
class FormulaSineWave : public Formula {
 public:
  //@{
  /** @name Arguments
    - amplitude: the amplitude (deviation from "average") of the wave
      (default: 1).
    - width: the width of the periodicity, w, in x units (default: 2\f$\pi\f$).
    - phase: shift the wave to the right in units of x (default: 0).
    - shift: shift the wave up in units of y (default: documented below).
   */
  explicit FormulaSineWave(argtype args = argtype());
  explicit FormulaSineWave(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Return the amplitude.
  double amplitude() const { return amplitude_; }

  /// Return the width.
  double width() const { return width_; }

  /// Return the phase
  double phase() const { return phase_; }

  // set the phase
  void set_phase(const double phase) { phase_ = phase; }

  /// Return the default value for the shift input argument.
  double default_shift() const { return 0.; }

  /// Return the shift
  double shift() const { return shift_; }

  // set the shift
  void set_shift(const double shift) { shift_ = shift; }

  /// Return the x-value of the minimum nearest to the point with the given
  /// x-value.
  double nearest_minimum(const double x) const;

  /// Same as above, but for maximum.
  double nearest_maximum(const double x) const;

  double evaluate(const double x) const override;
  double derivative(const double x) const override;
  double second_derivative(const double x) const;

  std::shared_ptr<Formula> create(std::istream& istr) const override {
    return std::make_shared<FormulaSineWave>(istr); }
  std::shared_ptr<Formula> create(argtype * args) const override {
    return std::make_shared<FormulaSineWave>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit FormulaSineWave(std::istream& istr);
  virtual ~FormulaSineWave() {}

  //@}
 private:
  double amplitude_, width_, phase_, shift_;
};

inline std::shared_ptr<FormulaSineWave> MakeFormulaSineWave(
    argtype args = argtype()) {
  return std::make_shared<FormulaSineWave>(args);
}

}  // namespace feasst

#endif  // FEASST_SHAPE_FORMULA_SINE_WAVE_H_
