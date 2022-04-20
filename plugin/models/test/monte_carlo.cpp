#include <memory>
#include "utils/test/utils.h"
#include "utils/include/io.h"
#include "math/include/accumulator.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/domain.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_transfer.h"
#include "monte_carlo/include/trial_rotate.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_add.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/metropolis.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model_intra.h"
#include "system/include/visit_model_cell.h"
#include "steppers/include/num_particles.h"
#include "steppers/include/log_and_movie.h"
//#include "models/include/lennard_jones_cut_shift.h"
#include "models/include/lennard_jones_force_shift.h"
#include "steppers/include/check_energy.h"

namespace feasst {

TEST(MonteCarlo, trimer) {
  MonteCarlo mc;
  mc.set(MakeRandomMT19937({{"seed", "123"}}));
  { { auto config = MakeConfiguration({{"cubic_box_length", "12"},
        {"particle_type", "../forcefield/trimer.fstprt"},
        {"add_particles_of_type0", "2"}});
      TrialSelectParticle sel;
      sel.select_particle(1, *config);
      const Position disp = Position().set_vector({4, 4, 4});
      config->displace_particle(sel.mobile(), disp);
      mc.add(config);
    }

    auto lj_wca = MakeLennardJonesForceShift();
    //auto lj_wca = MakeLennardJonesCutShift();
    ModelParams params = mc.system().configuration().model_params();
    lj_wca->set_wca(0, 1, &params);
    lj_wca->set_wca(1, 1, &params);
    //lj_wca->precompute(params);
    auto pot = MakePotential(lj_wca);
    pot->set(params);
    mc.add(pot);
  }
  mc.set(MakeThermoParams({{"beta", "4"}, {"chemical_potential", "-1"}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeTrialRotate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeTrialTransfer({{"particle_type", "0"}}));
  mc.add(MakeLogAndMovie({{"file_name", "tmp/trimer"}, {"trials_per", "1e2"}}));
  mc.add(MakeCheckEnergy({{"trials_per", "100"}, {"tolerance", "1e-10"}}));
  mc.attempt(1e3);
}

}  // namespace feasst
