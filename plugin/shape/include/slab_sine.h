
#ifndef FEASST_SHAPE_SLAB_SINE_H_
#define FEASST_SHAPE_SLAB_SINE_H_

#include <memory>
#include "utils/include/arguments.h"
#include "shape/include/formula_sine_wave.h"
#include "shape/include/shape_intersect.h"

namespace feasst {

/**
  Intersection of two parallel HalfSpaceSine.
  To impose 3-fold symmetry about the Cartesian origin and axes,
  The phases of the upper and lower HalfSpaceSine are adjusted.
  In particular, a quarter width of the sine wave is subtracted from the phase
  of the upper, and a quarter width is added to the phase of the lower.
  Thus, if FormulaSineWave has zero phase, then the crest of the upper and
  trough of the lower correspond with the zero in the wave dimension.
  Thus, the origin is the widest region when zero phase is given.
  If FormulaSineWave has half-width phase, then the origin is the most
  narrow region.
  Note that Domain::side_length in the wave_dimension must be an integer number
  of widths for periodicity.
 */
class SlabSine : public ShapeIntersect {
 public:
  /**
    args:
    - dimension: The slab surface is perpendicular to this dimensional axis.
    - wave_dimension : the wave travels along this dimension.
    - average_bound0: Set an average lower or upper value of the slab.
    - average_bound1: Set the second average bound, upper or lower, respectively.
    - FormulaSineWave arguments.
   */
  explicit SlabSine(argtype args = argtype());
  explicit SlabSine(argtype * args);

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Shape> create(std::istream& istr) const override {
    return std::make_shared<SlabSine>(istr); }
  std::shared_ptr<Shape> create(argtype * args) const override {
    return std::make_shared<SlabSine>(args); }
  explicit SlabSine(std::istream& istr);
  virtual ~SlabSine() {}
};

inline std::shared_ptr<SlabSine> MakeSlabSine(argtype args = argtype()) {
  return std::make_shared<SlabSine>(args);
}

}  // namespace feasst

#endif  // FEASST_SHAPE_SLAB_SINE_H_
