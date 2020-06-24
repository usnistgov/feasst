
#ifndef FEASST_SHAPE_SLAB_SINE_H_
#define FEASST_SHAPE_SLAB_SINE_H_

#include "utils/include/arguments.h"
#include "shape/include/shape.h"
#include "shape/include/formula_sine_wave.h"

namespace feasst {

/**
  Intersection of two parallel HalfSpaceSine.
  To impose 3-fold symmetry about the Cartesian origin and axes,
  The phases of the upper and lower HalfSpaceSine are adjusted.
  In particular, a quarter width of the sine wave is added to the phase of the
  upper, and a quarter width is subtracted from the phase of the lower.
  Thus, if FormulaSineWave has zero phase, then the trough of the upper and
  crest of the lower correspond with the origin.
  Thus, the origin is the most narrow region with zero phase.
  If FormulaSineWave has half-width phase, then the origin is the most
  expansive region.
  Note that Domain::side_length in the wave_dimension must be an integer number
  of widths for periodicity.
 */
class SlabSine : public Shape {
 public:
  /**
    args:
    - dimension: The slab surface is perpendicular to this dimensional axis.
    - wave_dimension : the wave travels along this dimension.
    - average_bound0: Set an average lower or upper value of the slab.
    - average_bound1: Set the second average bound, upper or lower, respectively.
   */
  SlabSine(std::shared_ptr<FormulaSineWave> sine_wave,
    const argtype &args = argtype());

  double nearest_distance(const Position& point) const override {
    return slab_->nearest_distance(point); }

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Shape> create(std::istream& istr) const override {
    return std::make_shared<SlabSine>(istr); }
  explicit SlabSine(std::istream& istr);
  virtual ~SlabSine() {}

 private:
  std::shared_ptr<Shape> slab_;
  Arguments args_;
};

inline std::shared_ptr<SlabSine> MakeSlabSine(
    std::shared_ptr<FormulaSineWave> sine_wave,
    const argtype& args = argtype()) {
  return std::make_shared<SlabSine>(sine_wave, args);
}

}  // namespace feasst

#endif  // FEASST_SHAPE_SLAB_SINE_H_
