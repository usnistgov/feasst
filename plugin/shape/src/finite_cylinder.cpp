#include "utils/include/io.h"
#include "utils/include/arguments.h"
#include "utils/include/debug.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "shape/include/half_space_tilted.h"
#include "shape/include/cylinder.h"
#include "shape/include/finite_cylinder.h"

namespace feasst {

FiniteCylinder::FiniteCylinder(argtype * args) : ShapeIntersect() {
  auto cyl = std::make_shared<Cylinder>(args);
  ASSERT(args->size() == 0, "unrecognized args: " << str(*args));
  auto first_endcap = MakeHalfSpaceTilted(cyl->first_point(), cyl->second_point());
  auto second_endcap = MakeHalfSpaceTilted(cyl->second_point(), cyl->first_point());
  set(MakeShapeIntersect(cyl, first_endcap), second_endcap);
}
FiniteCylinder::FiniteCylinder(argtype args) : FiniteCylinder(&args) {
  feasst_check_all_used(args);
}

}  // namespace feasst
