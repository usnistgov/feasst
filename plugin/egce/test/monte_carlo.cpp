#include "utils/test/utils.h"
#include "system/include/hard_sphere.h"
#include "monte_carlo/include/monte_carlo.h"
#include "ewald/test/system_example.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_rotate.h"
#include "steppers/include/movie.h"
#include "steppers/include/log.h"
#include "steppers/include/tuner.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/check_properties.h"
#include "steppers/include/cpu_time.h"
#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/macrostate_num_particles.h"
#include "flat_histogram/include/transition_matrix.h"
#include "steppers/include/criteria_writer.h"
#include "egce/include/a_equal_or_one_more_than_b.h"

namespace feasst {

MonteCarlo mc_rpm() {
  MonteCarlo mc;
  std::shared_ptr<Ewald> ewald;
  {
    System system;
    {
      Configuration config(MakeDomain({{"cubic_box_length", "12"}}), {
        {"particle_type0", "../plugin/ewald/forcefield/data.rpm_plus"},
        {"particle_type1", "../plugin/ewald/forcefield/data.rpm_minus"}
      });
      config.add_model_param("alpha", 4.8913043/config.domain()->min_side_length());
//      config.add_particle_of_type(0);
//      config.add_particle_of_type(1);
//      EXPECT_EQ(0, config.particle_type(0).type());
//      EXPECT_EQ(1, config.particle_type(1).type());
//      config.update_positions({{0, 0, 0},
//                               {2, 2, 2}});
//      INFO(config.particle(0).site(0).position().str());
//      INFO(config.particle(1).site(0).position().str());
      const int kmax_squared = 38;
      ewald = add_ewald_with(MakeHardSphere(), &system, kmax_squared);
      system.add(config);
    }
    mc.set(system);
  }
  //INFO(mc.system().configuration().particle_type(0).site(0).properties().str());
  //INFO(feasst_str(ewald->struct_fact_real_));
  const int steps_per = 1e2;
  mc.set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "-9"}}));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeTrialRotate({{"weight", "1."}, {"tunable_param", "1."}}));
  add_trial_transfer(&mc, {{"weight", "1."}, {"particle_type", "0"}});
  add_trial_transfer(&mc, {{"weight", "1."}, {"particle_type", "1"}});
  mc.add(MakeMovie({{"file_name", "tmp/rpm.xyz"}, {"steps_per", str(steps_per)}}));
  mc.add(MakeLog({{"file_name", "tmp/rpm_log.txt"}, {"steps_per", str(steps_per)}}));
  mc.add(MakeTuner({{"steps_per", str(steps_per)}}));
  mc.add(MakeCheckEnergy({{"tolerance", "1e-8"}, {"steps_per", str(steps_per)}}));
  mc.add(MakeCheckProperties({{"steps_per", str(steps_per)}}));
  mc.add(MakeCPUTime({{"steps_per", str(5*steps_per)}}));
  return mc;
}

TEST(MonteCarlo, rpm_egce) {
  MonteCarlo mc = mc_rpm();
  mc.set(MakeRandomMT19937({{"seed", "1574171557"}}));
  mc.attempt(1e3);
  test_serialize(mc);
  INFO(mc.analyze(2)->accumulator().str());
}

TEST(MonteCarlo, rpm_egce_fh) {
  MonteCarlo mc = mc_rpm();
  mc.set(MakeRandomMT19937({{"seed", "1574171557"}}));
  {
    auto criteria = MakeFlatHistogram({
      {"beta", "0.2"}, {"chemical_potential", "-250"}
    });
    {
      auto macro = MakeMacrostateNumParticles(
        Histogram({{"width", "1"}, {"max", "5"}}),
        {{"soft_max", "5"}}
      );
      macro->add(std::make_shared<AEqualOrOneMoreThanB>());
      criteria->set(macro);
    }
    criteria->set(MakeTransitionMatrix({
      {"min_sweeps", "10"},
      {"num_steps_to_update", "100"},
    }));
    mc.set(criteria);
  }
  mc.add(MakeCriteriaWriter({
    {"steps_per", str(1)},
    {"file_name", "tmp/rpmcrit.txt"},
  }));
  mc.attempt(1e3);
//  mc.run_until_complete();
  // test_serialize(mc);
}

}  // namespace feasst
