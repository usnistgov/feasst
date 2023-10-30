#include "utils/test/utils.h"
#include "utils/include/debug.h"
#include "utils/include/io.h"
#include "shape/include/half_space_tilted.h"
#include "shape/include/finite_cylinder.h"

namespace feasst {

TEST(FiniteCylinder, serialize) {
  auto finite_cylinder = MakeFiniteCylinder({{"radius", "5"},
    {"first_point", "f"}, {"f0", "0"}, {"f1", "0"}, {"f2", "0"},
    {"second_point", "s"}, {"s0", "0"}, {"s1", "0"}, {"s2", "2"}});
  std::shared_ptr<Shape> finite_cylinder2 = test_serialize<FiniteCylinder, Shape>(*finite_cylinder);

  EXPECT_TRUE (finite_cylinder2->is_inside(Position({0, 0, 1.49999}), 1.));
  EXPECT_FALSE(finite_cylinder2->is_inside(Position({0, 0, 1.50001}), 1.));
  EXPECT_TRUE (finite_cylinder2->is_inside(Position({0, 0, 0.50001}), 1.));
  EXPECT_FALSE(finite_cylinder2->is_inside(Position({0, 0, 0.49999}), 1.));
  EXPECT_TRUE (finite_cylinder2->is_inside(Position({0, 4.49999, 1}), 1.));
  EXPECT_FALSE(finite_cylinder2->is_inside(Position({0, 4.50001, 1}), 1.));

//  // visualize
//  std::vector<Position> grid = finite_cylinder2->grid(Position({6, 6, 4}),
//    Position({-6, -6, -4}));
//  for (const Position& pos : grid) {
//    std::cout << pos.coord(0) << " " << pos.coord(1) << " " << pos.coord(2) << std::endl;
//  }
}

}  // namespace feasst
