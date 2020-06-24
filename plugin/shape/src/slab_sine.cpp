#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "math/include/constants.h"
#include "shape/include/shape_intersect.h"
#include "shape/include/half_space_sine.h"
#include "shape/include/slab_sine.h"

namespace feasst {

class MapSlabSine {
 public:
  MapSlabSine() {
    auto obj = MakeSlabSine(MakeFormulaSineWave(), {
      {"dimension", "0"},
      {"wave_dimension", "1"},
      {"average_bound0", "1"},
      {"average_bound1", "-1"},
    });
    obj->deserialize_map()["SlabSine"] = obj;
  }
};

static MapSlabSine mapper_ = MapSlabSine();

SlabSine::SlabSine(std::shared_ptr<FormulaSineWave> sine_wave,
    const argtype &args) : Shape() {
  class_name_ = "SlabSine";
  args_.init(args);
  int dimension = args_.key("dimension").integer();
  int wave_dimension = args_.key("wave_dimension").integer();
  double upper = args_.key("average_bound0").dble();
  double lower = args_.key("average_bound1").dble();
  sort(&lower, &upper);
  ASSERT(upper - lower > NEAR_ZERO, "slab is infinitesimally thin");
  sine_wave->set_phase(sine_wave->phase() + 0.25*sine_wave->width());
  auto half0 = MakeHalfSpaceSine(sine_wave, {
    {"dimension", str(dimension)},
    {"wave_dimension", str(wave_dimension)},
    {"intersection", str(upper)},
    {"direction", str(-1)},
  });
  sine_wave->set_phase(sine_wave->phase() - 0.5*sine_wave->width());
  auto half1 = MakeHalfSpaceSine(sine_wave, {
    {"dimension", str(dimension)},
    {"wave_dimension", str(wave_dimension)},
    {"intersection", str(lower)},
    {"direction", str(1)},
  });
  slab_ = std::make_shared<ShapeIntersect>(half0, half1);
}

void SlabSine::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_shape_(ostr);
  feasst_serialize_version(1487, ostr);
  feasst_serialize(slab_, ostr);
}

SlabSine::SlabSine(std::istream& istr) : Shape(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(1487 == version, version);
  // feasst_deserialize(slab_, istr);
  // HWH for unknown reasons, the above doesn't work
  int existing;
  istr >> existing;
  if (existing != 0) {
    slab_ = slab_->deserialize(istr);
  }
}

}  // namespace feasst
