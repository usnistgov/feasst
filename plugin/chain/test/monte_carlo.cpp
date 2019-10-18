#include <memory>
#include "utils/test/utils.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_rotate.h"
#include "chain/include/trial_grow.h"
#include "chain/include/trial_pivot.h"
#include "chain/include/trial_crankshaft.h"
#include "chain/include/trial_reptate.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/test/monte_carlo_test.h"
#include "monte_carlo/include/criteria_metropolis.h"
#include "utils/include/utils_io.h"
#include "math/include/accumulator.h"
#include "system/test/system_test.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model_intra.h"
#include "system/include/visit_model_cell.h"
#include "steppers/include/log.h"
#include "steppers/include/movie.h"
#include "steppers/include/tuner.h"
#include "steppers/include/check.h"
#include "steppers/include/check_energy.h"
#include "chain/include/analyze_rigid_bonds.h"

namespace feasst {

Configuration config() {
  Configuration config({
    {"cubic_box_length", "12"},
    {"particle_type0", "../forcefield/data.chain10"},
    //{"particle_type0", "../plugin/chain/forcefield/data.chain50"},
    {"init_cells", "1."},
  });
  config.add_particle_of_type(0);
  return config;
}

Potential lj() {
  Potential lj;
  lj.set_model(std::make_shared<ModelLJ>());
  return lj;
}

Potential lj_dual_cut(const Configuration config) {
  Potential lj_dual_cut = lj();
  lj_dual_cut.set_model_params(config);
  lj_dual_cut.set_model_param("cutoff", 0, 1);
  lj_dual_cut.set_visit_model(std::make_shared<VisitModelCell>());
  return lj_dual_cut;
}

Potential lj_intra() {
  Potential lj_intra = lj();
  auto visitor = std::make_shared<VisitModelIntra>();
  visitor->set_intra_cut(1);
  lj_intra.set_visit_model(visitor);
  return lj_intra;
}

Potential lj_intra_dual_cut(const Configuration config) {
  Potential lj_intra_dual_cut = lj_intra();
  lj_intra_dual_cut.set_model_params(config);
  lj_intra_dual_cut.set_model_param("cutoff", 0, 1);
  return lj_intra_dual_cut;
}

Potential lrc() {
  Potential lrc;
  lrc.set_visit_model(std::make_shared<LongRangeCorrections>());
  return lrc;
}

System chain_system() {
  System system;
  system.add(config());
  system.add_to_unoptimized(lj());
  system.add_to_reference(lj_dual_cut(system.configuration()));
  system.add_to_unoptimized(lj_intra());
  system.add_to_reference(lj_intra_dual_cut(system.configuration()));
  system.add_to_unoptimized(lrc());
  return system;
}

TEST(MonteCarlo, chain) {
  MonteCarlo mc;
  mc.set(chain_system());
  mc.set(MakeCriteriaMetropolis({{"beta", "1"}, {"chemical_potential", "1."}}));
  mc.seek_num_particles(1);
  mc.add(MakeTrialTranslate({
    {"weight", "1."},
    {"tunable_param", "1."},
    {"reference_index", "0"},
    {"num_steps", "2"},
  }));
  mc.add(MakeTrialRotate({
    {"weight", "1."},
    {"tunable_param", "20."},
    {"reference_index", "0"},
    {"num_steps", "2"},
  }));
  mc.add(MakeTrialPivot({
    {"weight", "1."},
    {"tunable_param", "20."},
    {"max_length", "30"},
    {"reference_index", "0"},
    {"num_steps", "2"},
  }));
  mc.add(MakeTrialReptate({
    {"weight", "1."},
//    {"reference_index", "0"},
//    {"num_steps", "2"},
  }));
  mc.add(MakeTrialCrankshaft({
    {"weight", "1."},
    {"tunable_param", "25."},
    {"max_length", "5."},
    {"reference_index", "0"},
    {"num_steps", "2"},
  }));
  mc.add(MakeTrialGrowLinear(
    std::make_shared<TrialComputeMove>(),
    {
//      {"weight", "0.1"},
      {"particle_type", "0"},
      {"num_steps", "3"},
      {"reference", "0"},
    }
  ));
  const int steps_per = 1e2;
  mc.add(MakeLog({
    {"steps_per", str(steps_per)},
    {"file_name", "tmp/chainlog.txt"},
  }));
  mc.add(MakeMovie({
    {"steps_per", str(steps_per)},
    {"file_name", "tmp/chain10movie.xyz"},
  }));
  mc.add(MakeCheckEnergy({
    {"steps_per", str(steps_per)},
    {"tolerance", "1e-10"},
  }));
  mc.add(MakeTuner({{"steps_per", str(steps_per)}}));
  mc.add(MakeAnalyzeRigidBonds({{"steps_per", str(steps_per)}}));
  mc.attempt(3e2);

  MonteCarlo mc2 = test_serialize(mc);
  EXPECT_EQ(mc2.analyzers().size(), 3);
}

}  // namespace feasst
