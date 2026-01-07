
#include "utils/test/utils.h"
#include "math/include/recursive_table.h"

namespace feasst {

TEST(RecursiveTable1D, serialize) {
  auto obj = std::make_unique<RecursiveTable1D>(argtype({{"num", "2"}}));
  auto obj2 = test_serialize_unique(*obj);
}

/**
  |                |   o...o   o   |               |
 z=0             z=1/3           z=2/3            z=1
f(z)=1         f(z)=0.5         f(z)=0           f(z)=-0.1

      second level -> .7  .3  .8
  third level between 0.7 and 0.3-> 1.1 -0.5
 */
TEST(RecursiveTable1D, linear) {
  auto tab = std::make_unique<RecursiveTable1D>(argtype({{"num", "4"}}));
  tab->set_data(0, 1.);
  tab->set_data(1, 0.5);
  tab->set_data(2, 0.0);
  tab->set_data(3, -0.1);
  auto tab2 = std::make_unique<RecursiveTable1D>(argtype({{"num", "5"}}));
  tab2->set_data(0, 0.5);
  tab2->set_data(1, 0.7);
  tab2->set_data(2, 0.3);
  tab2->set_data(3, 0.8);
  tab2->set_data(4, 0.0);
  auto tab3 = std::make_unique<RecursiveTable1D>(argtype({{"num", "4"}}));
  tab3->set_data(0, 0.7);
  tab3->set_data(1, 1.1);
  tab3->set_data(2, -0.5);
  tab3->set_data(3, 0.3);
  tab2->insert(1, *tab3);
  tab->insert(1, *tab2);
  EXPECT_NEAR(tab->linear_interpolation(2./3. + 1./6.), -0.05, 1e-12);
  EXPECT_NEAR(tab->linear_interpolation(0.4), 0.5+0.2/0.25*0.2, 1e-12);
  EXPECT_NEAR(tab->linear_interpolation(0.45), 0.78, 1e-12);
}

}  // namespace feasst
