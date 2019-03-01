#include <memory>
#include <gtest/gtest.h>
#include "core/include/trial_translate.h"
#include "core/include/trial_rotate.h"
#include "chain/include/trial_pivot.h"
#include "chain/include/trial_crankshaft.h"
#include "chain/include/trial_regrow.h"
#include "core/include/trial_transfer.h"
#include "core/include/monte_carlo.h"
#include "core/include/criteria_metropolis.h"
#include "core/include/utils_io.h"
#include "core/include/accumulator.h"
#include "core/test/system_test.h"
#include "core/include/long_range_corrections.h"
#include "core/include/visit_model_intra.h"
#include "core/include/visit_model_cell.h"
#include "chain/include/analyze_rigid_bonds.h"

namespace feasst {

TEST(MonteCarlo, chain) {
  seed_random_by_date();
  seed_random(1551407871);
  MonteCarlo mc;

  { System system;
    { Configuration config;
      config.set_domain(Domain().set_cubic(12));
      config.add_particle_type("../forcefield/data.chain10");
      system.add(config); }

    { Potential potential;
      potential.set_model(std::make_shared<ModelLJ>());
      potential.set_visit_model(std::make_shared<VisitModel>());
      potential.set_model_params(system.configuration());
      system.add_to_unoptimized(potential); }

    { Potential potential;
      potential.set_model(std::make_shared<ModelHardSphere>());
      auto visitor = std::make_shared<VisitModelIntra>();
      visitor->set_intra_cut(1);
      potential.set_visit_model(visitor);
      potential.set_model_params(system.configuration());
      potential.set_model_param("cutoff", 0, 1.);
      system.add_to_unoptimized(potential); }

    { Potential lrc;
      lrc.set_visit_model(std::make_shared<LongRangeCorrections>());
      lrc.set_model_params(system.configuration());
      system.add_to_unoptimized(lrc); }

    mc.set(system);
  }

  mc.set(MakeCriteriaMetropolis({{"beta", "1.2"}, {"add_activity", "1."}}));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"max_move", "1."}}));

  { auto trial = std::make_shared<TrialRotate>();
    trial->set_weight(1.);
    trial->set_max_move(90.);
    mc.add(trial); }

  { auto trial = std::make_shared<TrialPivot>();
    trial->set_weight(1.);
    trial->set_max_move(90.);
    mc.add(trial); }

  { auto trial = std::make_shared<TrialCrankshaft>();
    trial->set_weight(1.);
    trial->set_max_move(90.);
    trial->set_tunable_percent_change(0.1);
    mc.add(trial); }

  { auto trial = std::make_shared<TrialRegrow>();
    trial->set_weight(1.);
    mc.add(trial); }

  mc.seek_num_particles(1);
  const int num_periodic = 1e3;

  mc.add(MakeLog(
   {{"steps_per", str(num_periodic)},
    {"file_name", "tmp/chainlog.txt"}}));
  mc.add(MakeMovie(
   {{"steps_per", str(num_periodic)},
    {"file_name", "tmp/chain10movie.xyz"}}));
  mc.add(MakeEnergyCheck(
   {{"steps_per", str(num_periodic)},
    {"tolerance", "1e-10"}}));
  mc.add(MakeTuner({{"steps_per", str(num_periodic)}}));
  mc.add(MakeAnalyzeRigidBonds({{"steps_per", str(num_periodic)}}));

  mc.attempt(1e4);
}

}  // namespace feasst
