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
#include "core/include/criteria_flat_histogram.h"
#include "core/include/macrostate_num_particles.h"
#include "core/include/bias_wang_landau.h"
#include "core/include/histogram.h"
#include "core/include/utils_io.h"
#include "core/include/accumulator.h"
#include "core/test/system_test.h"
#include "core/include/long_range_corrections.h"
#include "core/include/visit_model_intra.h"
#include "core/include/visit_model_cell.h"

namespace feasst {

TEST(MonteCarlo, chain) {
  seed_random_by_date();
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

  { auto criteria = std::make_shared<CriteriaMetropolis>();
    criteria->set_beta(1.2);
    criteria->add_activity(1.);
    mc.set(criteria); }

  { auto trial = std::make_shared<TrialTranslate>();
    trial->set_weight(1.);
    trial->set_max_move(1.);
    trial->set_max_move_bounds(mc.system().configuration().domain());
    mc.add(trial); }

  { auto trial = std::make_shared<TrialRotate>();
    trial->set_weight(1.);
    trial->set_max_move(90.);
    mc.add(trial); }

  { auto trial = std::make_shared<TrialPivot>();
    trial->set_weight(1.);
    trial->set_max_move(90.);
    mc.add(trial); }

  { auto trial = std::make_shared<TrialCrankShaft>();
    trial->set_weight(1.);
    trial->set_max_move(90.);
    trial->set_tunable_percent_change(0.1);
    mc.add(trial); }

  { auto trial = std::make_shared<TrialRegrow>();
    trial->set_weight(1.);
    mc.add(trial); }

  mc.seek_num_particles(1);
  const int num_periodic = 1e3;

  { auto log = std::make_shared<Log>();
    log->set_steps_per_write(num_periodic);
    log->set_file_name("tmp/log.txt");
    mc.add(log); }

  { auto movie = std::make_shared<Movie>();
    movie->set_steps_per(num_periodic);
    movie->set_file_name("tmp/chain10movie.xyz");
    mc.add(movie); }

  { auto checker = std::make_shared<EnergyCheck>();
    checker->set_steps_per_update(num_periodic);
    checker->set_tolerance(1e-10);
    mc.add(checker); }

  { auto tuner = std::make_shared<Tuner>();
    tuner->set_steps_per_update(num_periodic);
    mc.add(tuner); }

  { auto bondcheck = std::make_shared<RigidBondChecker>();
    bondcheck->set_steps_per(num_periodic);
    mc.add(bondcheck); }

  mc.attempt(1e4);
}

}  // namespace feasst
