#include <cmath>
#include "utils/include/utils_io.h"
#include "utils/include/debug.h"
#include "system/include/lennard_jones.h"
#include "system/include/hard_sphere.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model_bond.h"
#include "system/include/model_two_body_factory.h"
#include "monte_carlo/include/trial_rotate.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/metropolis.h"
#include "steppers/include/movie.h"
#include "steppers/include/log.h"
#include "steppers/include/tuner.h"
#include "steppers/include/check_energy.h"
#include "ewald/include/utils.h"
#include "ewald/include/ewald.h"
#include "ewald/include/charge_screened.h"
#include "ewald/include/charge_screened_intra.h"
#include "ewald/include/charge_self.h"
#include "ewald/include/trial_add_multiple.h"
#include "ewald/include/trial_remove_multiple.h"

namespace feasst {

void add_trial_transfer_multiple(MonteCarlo * mc, const argtype& args) {
  mc->add(MakeTrialAddMultiple(args));
  mc->add(MakeTrialRemoveMultiple(args));
}

void spce(MonteCarlo * mc, const argtype& args) {
  Arguments args_(args);
  { Configuration config(
      MakeDomain({{"cubic_box_length",
                   args_.key("cubic_box_length").dflt("24.8586887").str()}}),
      {{"particle_type", install_dir() + "/" + args_.key("particle").dflt("forcefield/data.spce").str()},
       {"physical_constants", args_.key("physical_constants").dflt("CODATA2018").str()}}
    );
    if (args_.key("xyz_file").used()) {
      FileXYZ().load(args_.str(), &config);
    }
    mc->add(config);
  }
  mc->add(Potential(
    MakeEwald({{"kmax_squared", args_.key("kmax_squared").dflt("38").str()},
               {"alpha",
      str(args_.key("alphaL").dflt("5.6").dble()/
          mc->configuration().domain().min_side_length())}})));
  mc->add(Potential(MakeModelTwoBodyFactory({MakeLennardJones(),
                                             MakeChargeScreened()})));
  mc->add(Potential(MakeChargeScreenedIntra(), MakeVisitModelBond()));
  mc->add(Potential(MakeChargeSelf()));
  mc->add(Potential(MakeLongRangeCorrections()));
  const double R = mc->configuration().physical_constants().ideal_gas_constant();
  const double temp = args_.key("temperature").dflt("525").dble()*R/1000;  // KJ/mol
  const std::string mu =
    args_.key("chemical_potential").dflt("-35.294567543492").str();
  mc->set(MakeMetropolis({{"beta", str(1/temp)},
                          {"chemical_potential", mu}}));
  mc->add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "0.275"}}));
  mc->add(MakeTrialRotate({{"weight", "1."}, {"tunable_param", "50."}}));
}

void rpm(MonteCarlo * mc, const argtype& args) {
  Arguments args_(args);
  { Configuration config(
      MakeDomain({{"cubic_box_length",
                   args_.key("cubic_box_length").dflt("12").str()}}),
      {{"particle_type0", install_dir() + "/" + args_.key("particle0").dflt("plugin/ewald/forcefield/data.rpm_plus").str()},
       {"particle_type1", install_dir() + "/" + args_.key("particle1").dflt("plugin/ewald/forcefield/data.rpm_minus").str()}}
    );
    if (args_.key("cutoff").used()) {
      const double cutoff = args_.dble();
      config.set_model_param("cutoff", 0, cutoff);
      config.set_model_param("cutoff", 1, cutoff);
    }
    const double charge_conversion = config.physical_constants().charge_conversion();
    const Charge& charge = config.model_params().charge();
    for (int type = 0; type < config.num_particle_types(); ++type) {
      double cval = charge.value(type);
      if (type == 0) {
        cval *= args_.key("charge_ratio").dflt("1").dble();
      }
      config.set_model_param("charge", type,
        cval/std::sqrt(charge_conversion));
    }
    if (args_.key("delta").used()) {
      const double delta = args_.dble();
      const Sigma& sigma = config.model_params().sigma();
      config.set_model_param("sigma", 0, sigma.value(0) + delta);
      config.set_model_param("sigma", 1, sigma.value(1) - delta);
    }
    mc->add(config);
  }
  mc->add(Potential(
    MakeEwald({{"kmax_squared", args_.key("kmax_squared").dflt("38").str()},
               {"alpha",
      str(args_.key("alphaL").dflt("5.6").dble()/
          mc->configuration().domain().min_side_length())}})));
  mc->add(Potential(MakeModelTwoBodyFactory({MakeHardSphere(),
                                             MakeChargeScreened()})));
  mc->add(Potential(MakeChargeSelf()));
  const double temp = args_.key("temperature").dflt("0.047899460618081").dble();
  mc->set(MakeMetropolis({{"beta", str(1./temp)},
    {"chemical_potential", str(args_.key("beta_mu").dflt("-13.94").dble()/temp)}}));
  mc->add(MakeTrialTranslate({{"weight", "0.25"}, {"tunable_param", "0.1"}}));
}

}  // namespace feasst
