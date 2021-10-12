#include "utils/test/utils.h"
#include "system/test/sys_utils.h"
#include "monte_carlo/test/perturb_test.h"
#include "cluster/include/energy_map_all.h"

namespace feasst {

TEST(Perturb, RevertEnergyMap) {
  System system = two_particle_system();
  system.set_unoptimized(0,
    MakePotential(MakeLennardJones(),
                  MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
  test_revert(&system);
}

}  // namespace feasst
