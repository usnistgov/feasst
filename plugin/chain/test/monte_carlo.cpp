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
  seed_random();
  MonteCarlo mc;

  { System system;
    { Configuration config({
        {"cubic_box_length", "12"},
        {"particle_type0", "../forcefield/data.chain10"},
        {"init_cells", "1."},
      });
      config.add_particle_of_type(0);
      system.add(config); }

    { Potential potential;
      potential.set_model(std::make_shared<ModelLJ>());
      potential.set_visit_model(std::make_shared<VisitModel>());
      system.add_to_unoptimized(potential);
      potential.set_model_params(system.configuration());
      potential.set_model_param("cutoff", 0, 1);
      potential.set_visit_model(std::make_shared<VisitModelCell>());
      system.add_to_reference(potential);
    }

    { Potential potential;
      potential.set_model(std::make_shared<ModelHardSphere>());
      auto visitor = std::make_shared<VisitModelIntra>();
      visitor->set_intra_cut(1);
      potential.set_visit_model(visitor);
      potential.set_model_params(system.configuration());
      potential.set_model_param("cutoff", 0, 1.);
      system.add_to_unoptimized(potential);
    }

    { Potential lrc;
      lrc.set_visit_model(std::make_shared<LongRangeCorrections>());
      lrc.set_model_params(system.configuration());
      system.add_to_unoptimized(lrc);
    }

    mc.set(system);
  }

  mc.set(MakeCriteriaMetropolis({{"beta", "1.2"}, {"add_activity", "1."}}));
  mc.seek_num_particles(2);
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"max_move", "1."}}));
  mc.add(MakeTrialRotate({{"weight", "1."}, {"max_move", "20."}}));
  mc.add(MakeTrialPivot({{"weight", "1."}, {"max_move", "20."}}));
  mc.add(MakeTrialCrankshaft({{"weight", "1."}, {"max_move", "1."}}));
  mc.add(MakeTrialRegrow({
    {"weight", "1."},
    {"num_steps", "5"},
    {"reference", "0"},
  }));
  const int steps_per = 1e2;
  mc.add(MakeLog({
    {"steps_per", str(steps_per)},
    {"file_name", "tmp/chainlog.txt"},
  }));
  mc.add(MakeMovie({
    {"steps_per", str(steps_per)},
    {"file_name", "tmp/chain10movie.xyz"},
  }));
  mc.add(MakeEnergyCheck({
    {"steps_per", str(steps_per)},
    {"tolerance", "1e-10"},
  }));
  mc.add(MakeTuner({{"steps_per", str(steps_per)}}));
  mc.add(MakeAnalyzeRigidBonds({{"steps_per", str(steps_per)}}));
  mc.attempt(1e3);
}

}  // namespace feasst
