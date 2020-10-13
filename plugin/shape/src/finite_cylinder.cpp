#include "utils/include/io.h"
#include "utils/include/debug.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "shape/include/half_space_tilted.h"
#include "shape/include/cylinder.h"
#include "shape/include/finite_cylinder.h"

namespace feasst {

FiniteCylinder::FiniteCylinder(const argtype &args,
    const Position& point0,
    const Position& point1) : ShapeIntersect() {
  set(MakeShapeIntersect(MakeCylinder(args, point0, point1),
                         MakeHalfSpaceTilted(point0, point1)),
    MakeHalfSpaceTilted(point1, point0));
}

}  // namespace feasst
