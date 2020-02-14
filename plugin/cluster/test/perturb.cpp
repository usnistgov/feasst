#include <gtest/gtest.h>
#include "monte_carlo/test/perturb_test.h"
#include "cluster/include/energy_map_all.h"

namespace feasst {

TEST(Perturb, RevertEnergyMap) {
  System system = default_system();
  system.set_unoptimized(0,
    Potential(MakeLennardJones(),
              MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
  test_revert(&system);
}

}  // namespace feasst
