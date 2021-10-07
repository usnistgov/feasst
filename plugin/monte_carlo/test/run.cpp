
#include "utils/test/utils.h"
#include "system/include/lennard_jones.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/run.h"

namespace feasst {

TEST(AddReference, lj) {
  MonteCarlo mc;
  mc.add(MakeConfiguration({{"cubic_box_length", "8"}, {"particle_type0", "../forcefield/lj.fstprt"}}));
  mc.add(MakePotential(MakeLennardJones()));
  mc.run(MakeAddReference({{"cutoff", "1"}, {"use_cell", "true"}}));
  EXPECT_EQ(mc.system().num_references(), 1);
  EXPECT_EQ(mc.system().reference(0, 0).visit_model().class_name(), "VisitModelCell");
}

}  // namespace feasst
