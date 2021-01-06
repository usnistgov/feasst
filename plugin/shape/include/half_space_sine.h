
#ifndef FEASST_SHAPE_HALF_SPACE_SINE_H_
#define FEASST_SHAPE_HALF_SPACE_SINE_H_

#include "utils/include/arguments.h"
#include "shape/include/formula_sine_wave.h"
#include "shape/include/half_space.h"

namespace feasst {

class Solver;
class SineDistDeriv;

/**
  Similar to HalfSpace, which divides space by a plane (or line in 2D),
  except that space is divided by a sine wave instead.
  FormulaSineWave describes the surface.
  Note that the "y" coordinate refers to the one specified in the
  "dimension" input argument.
  Thus, the "intersect" argument in HalfSpace overrides the "shift" argument in
  FormulaSineWave.
  The "x" coordinate is given by wave_dimension.

  To find the nearest distance of the sine wall from the point,
  the squared distance, \f$D^2\f$ of the point (e,f) to the sine wave,

  \f$D^2 = (x-e)^2 + (y-f)^2\f$

  should be a minimum. Thus, the derivative should be zero.

  \f$\frac{d D^2}{2dx} = 0 = x-e+(y-f)*\frac{dy}{dx}

  For this implementation, we begin by finding the limits of the nearest
  half-wave of interest, in order to bracket the possible solutions to this
  non-linear equation.
  When f > y(e), the "x"-values of the half-wave of interest is given by the
  nearest minimum to e, and then the nearest maximum on the other side of e.
  Otherwise, when f < y(e), switch the minimum and maximum in above.
  When f = y(e), (e,f) is on the curve and your nearest distance is 0.
 */
class HalfSpaceSine : public HalfSpace {
 public:
  /**
    args:
    - wave_dimension : the wave travels along this dimension.
   */
  HalfSpaceSine(std::shared_ptr<FormulaSineWave> sine_wave,
    const argtype &args = argtype());

  /// Return the sine wave formula.
  const FormulaSineWave& sine_wave() const { return sine_wave_; }

  /// Return the dimension along which the wave travels.
  int wave_dimension() const { return wave_dimension_; }

  bool is_inside(const Position& point) const override;
  bool is_inside(const Position& point, const double diameter) const override;
  double nearest_distance(const Position& point) const override;

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Shape> create(std::istream& istr) const override {
    return std::make_shared<HalfSpaceSine>(istr); }
  explicit HalfSpaceSine(std::istream& istr);
  virtual ~HalfSpaceSine() {}

 private:
  int wave_dimension_;
  FormulaSineWave sine_wave_;
  std::shared_ptr<Solver> solver_;

  //temporary
  std::shared_ptr<SineDistDeriv> sine_dist_deriv_;
  void init_sine_dist_deriv_();
};

inline std::shared_ptr<HalfSpaceSine> MakeHalfSpaceSine(
    std::shared_ptr<FormulaSineWave> sine_wave,
    const argtype& args = argtype()) {
  return std::make_shared<HalfSpaceSine>(sine_wave, args);
}

}  // namespace feasst

#endif  // FEASST_SHAPE_HALF_SPACE_SINE_H_
