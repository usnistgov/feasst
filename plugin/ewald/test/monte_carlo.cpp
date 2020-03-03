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
#include "steppers/include/check_properties.h"
#include "steppers/include/cpu_time.h"

namespace feasst {

TEST(MonteCarlo, spce) {
  MonteCarlo mc;
  mc.set(MakeRandomMT19937({{"seed", "time"}}));
  // mc.set(MakeRandomMT19937({{"seed", "1234"}}));
  // mc.set(MakeRandomMT19937({{"seed", "1572272377"}}));
  // mc.set(MakeRandomMT19937({{"seed", "1574171557"}}));
  std::shared_ptr<Ewald> ewald;
  {
    System system;
    {
      Configuration config(MakeDomain({{"cubic_box_length", "24.8586887"}}),
        {{"particle_type", "../forcefield/data.spce"}});
      system.add(config);
    }
//    ewald = add_ewald_with(MakeLennardJones(), &system, kmax_squared);
    system.add(Potential(
      MakeEwald({{"kmax_squared", "3"},
                 {"alpha", str(5.6/system.configuration().domain()->min_side_length())}}),
                       {{"prevent_cache", "true"}}));
    system.add(Potential(MakeModelTwoBodyFactory({MakeLennardJones(),
                                               MakeChargeScreened()})));
    //system.add(Potential(MakeChargeScreenedIntra(),
    //                   MakeVisitModelIntra({{"cutoff", "0"}})));
    system.add(Potential(MakeChargeScreenedIntra(), MakeVisitModelBond()));
    system.add(Potential(MakeChargeSelf()));
    system.add(Potential(MakeLongRangeCorrections()));
    mc.set(system);
  }
  //INFO(mc.system().configuration().particle_type(0).site(0).properties().str());
  //INFO(feasst_str(ewald->struct_fact_real_));
  const int steps_per = 1e2;
  mc.set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "-9"}}));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeTrialRotate({{"weight", "1."}, {"tunable_param", "1."}}));
  // mc.add(MakeTrialRemove({{"weight", "1."}, {"particle_type", "0"}}));
  add_trial_transfer(&mc, {{"weight", "1."}, {"particle_type", "0"}});
  mc.add(MakeMovie({{"file_name", "tmp/spce.xyz"}, {"steps_per", str(steps_per)}}));
  mc.add(MakeLog({{"file_name", "tmp/spce_log.txt"}, {"steps_per", str(steps_per)}}));
  mc.add(MakeTuner({{"steps_per", str(1e2)}}));
  mc.add(MakeCheckEnergy({{"tolerance", "1e-8"}, {"steps_per", str(steps_per)}}));
  mc.add(MakeCheckProperties({{"steps_per", str(steps_per)}}));
  mc.add(MakeCPUTime({{"steps_per", str(5*steps_per)}}));

  //mc.seek_num_particles(2);
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
  INFO(mc.analyze(2)->accumulator().str());
}

}  // namespace feasst
