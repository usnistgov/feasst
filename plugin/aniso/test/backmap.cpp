#include "utils/test/utils.h"
#include "configuration/include/domain.h"
#include "configuration/include/physical_constants.h"
#include "monte_carlo/include/monte_carlo.h"
#include "aniso/include/tabulate_two_rigid_body_3D.h"
#include "aniso/include/backmap.h"

namespace feasst {

TEST(MonteCarlo, backmap_LONG) {
  WARN("update for more simple backmap with 2d L");
  auto mc = MakeMonteCarlo({{
    {"Configuration", {
      {"particle_type0", "../plugin/steppers/test/data/mab.txt"},
      {"xyz_euler_file", "../plugin/steppers/test/data/nvt0.xyze"},
    }},
    {"Potential", {{"VisitModel", "DontVisitModel"}}},
    {"ThermoParams", {{"beta", "1"}, {"chemical_potential", "1"}}},
    {"Metropolis", {{}}},
    //{"Metropolis", {{"trials_per_cycle", "10"}, {"cycles_to_complete", "10"}}},
    {"ReadConfigFromFile", {{"input_file", "../plugin/steppers/test/data/nvt0.xyze"}, {"euler", "true"}}},
    {"Backmap", {{"trials_per_write", "1"}, {"output_file", "tmp/backmap.xyz"}, {"site0", "0"}, {"fstprt0", "../plugin/aniso/test/data/fc.txt"}, {"site1", "3"}, {"fstprt1", "../plugin/aniso/test/data/fc.txt"}}},
    {"Run", {{"until", "complete"}}},
  }}, true);
}

}  // namespace feasst
