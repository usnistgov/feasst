#include "utils/test/utils.h"
#include "monte_carlo/include/monte_carlo.h"
#include "ewald/test/system_example.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_rotate.h"
#include "steppers/include/movie.h"
#include "steppers/include/log.h"
#include "steppers/include/tuner.h"
#include "steppers/include/check_energy.h"

namespace feasst {

TEST(MonteCarlo, spce) {
  MonteCarlo mc;
  // mc.set(MakeRandomMT19937({{"seed", "1234"}}));
  // mc.set(MakeRandomMT19937({{"seed", "1572272377"}}));
  std::shared_ptr<Ewald> ewald;
  {
    System system;
    {
      Configuration config({
        {"cubic_box_length", "24.8586887"},
        {"particle_type", "../forcefield/data.spce"}
      });
      config.add_model_param("alpha", 5.6/config.domain().min_side_length());
//      config.add_particle_of_type(0);
      const int kmax_squared = 3;
      ewald = add_ewald_with(MakeModelLJ(), &config, &system, kmax_squared);
      system.add(config);
    }
    mc.set(system);
  }
  //INFO(mc.system().configuration().particle_type(0).site(0).properties().str());
  //INFO(feasst_str(ewald->struct_fact_real_));
  const int steps_per = 1e1;
  mc.set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "-9"}}));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeTrialRotate({{"weight", "1."}, {"tunable_param", "1."}}));
  // mc.add(MakeTrialRemove({{"weight", "1."}, {"particle_type", "0"}}));
  add_trial_transfer(&mc, {{"weight", "1."}, {"particle_type", "0"}});
  mc.add(MakeMovie({{"file_name", "tmp/spce.xyz"}, {"steps_per", str(steps_per)}}));
  mc.add(MakeLog({{"file_name", "tmp/spce_log.txt"}, {"steps_per", str(steps_per)}}));
  mc.add(MakeTuner({{"steps_per", str(steps_per)}}));
  mc.add(MakeCheckEnergy({{"tolerance", "1e-8"}, {"steps_per", str(steps_per)}}));

  // Theres something wrong with MC and Ewald
  mc.seek_num_particles(2);
  // INFO(mc.system().configuration().num_particles());
//  INFO("props " << mc.system().configuration().particle(0).site(1).properties().str());
  mc.attempt(1e3);

//  INFO(ewald->num_vectors());
//  INFO(feasst_str(ewald->struct_fact_real_));
//  std::vector<double> reals = ewald->struct_fact_real_;
//  std::vector<double> imags = ewald->struct_fact_imag_;
//  INFO("props " << mc.system().configuration().particle(0).site(1).properties().str());
//  mc.get_system()->energy();
//  INFO(feasst_str(ewald->struct_fact_real_));
////  EXPECT_EQ(reals, ewald->struct_fact_real_);
//  EXPECT_TRUE(is_equal(reals, ewald->struct_fact_real_));
//  EXPECT_TRUE(is_equal(imags, ewald->struct_fact_imag_));
//  INFO("props " << mc.system().configuration().particle(0).site(1).properties().str());

  test_serialize(mc);
}

}  // namespace feasst
