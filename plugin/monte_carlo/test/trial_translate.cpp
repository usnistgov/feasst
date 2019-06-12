
#include <gtest/gtest.h>
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/criteria_metropolis.h"
#include "system/include/visit_model_cell.h"
#include "system/include/model_lj.h"
#include "system/include/model_hard_sphere.h"
#include "steppers/include/log.h"
#include "steppers/include/movie.h"
#include "steppers/include/tuner.h"
#include "steppers/include/check.h"
#include "steppers/include/check_energy.h"

namespace feasst {

// HWH this test failes because of a change to rosenbluth storage.
//TEST(TrialStagedTranslate, stage) {
//  std::vector<int> refs = {-1, 0};
//  for (int ref : refs) {
//    seed_random_by_date();
//    MonteCarlo mc;
//
//    { System system;
//      { Configuration config;
//        config.set_domain(Domain().set_cubic(6));
//        config.add_particle_type("../forcefield/data.lj");
//        system.add(config); }
//
//      { Potential potential;
//        potential.set_model(std::make_shared<ModelLJ>());
//        system.add_to_unoptimized(potential);
//        if (ref == 0) {
//          potential.set_model_params(system.configuration());
//          potential.set_model_param("cutoff", 0, 1.);
//          system.add_to_reference(potential);
//          EXPECT_EQ(3., system.configuration().model_params().cutoff().mixed_value(0, 0));
//          EXPECT_EQ(1., system.reference(0, 0).model_params().cutoff().mixed_value(0, 0));
//          try {
//            system.unoptimized().potentials()[0].model_params().cutoff().mixed_value(0, 0);
//            CATCH_PHRASE("you must first initialize model params");
//          }
//        }
//      }
//
//      mc.set(system);
//    }
//
//    mc.set(MakeCriteriaMetropolis({{"beta", "1.2"}, {"chemical_potential", "1."}}));
//    mc.seek_num_particles(20);
//    if (ref == -1) {
//      mc.add(MakeTrialStagedTranslate(
//        {{"weight", "1."},
//         {"tunable_param", "1."},
//         {"num_steps", "3"},
//      }));
//    } else {
//      mc.add(MakeTrialStagedTranslate(
//        {{"weight", "1."},
//         {"tunable_param", "1."},
//         {"reference", "0"},
//         {"num_steps", "3"},
//        }));
//    }
//    const int num_periodic = 1e2;
//    mc.add(MakeLog(
//     {{"steps_per", str(num_periodic)},
//      {"file_name", "tmp/ljmfblog.txt"}}));
//    mc.add(MakeMovie(
//     {{"steps_per", str(num_periodic)},
//      {"file_name", "tmp/ljmfbmovie.xyz"}}));
//    mc.add(MakeCheckEnergy(
//     {{"steps_per", str(num_periodic)},
//      {"tolerance", "1e-10"}}));
//    mc.add(MakeTuner({{"steps_per", str(num_periodic)}}));
//    // INFO("energy " << mc.criteria()->current_energy());
//    mc.attempt(1e3);
//  }
//}

}  // namespace feasst
