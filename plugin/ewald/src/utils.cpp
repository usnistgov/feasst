#include <cmath>
#include "utils/include/io.h"
#include "utils/include/debug.h"
#include "math/include/constants.h"
#include "configuration/include/domain.h"
#include "system/include/lennard_jones.h"
#include "system/include/hard_sphere.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model_bond.h"
#include "system/include/model_two_body_factory.h"
#include "system/include/visit_model_cell.h"
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

double kelvin2kJpermol(const double kelvin) {
  ModelParams model_params;
  const double R = model_params.physical_constants().ideal_gas_constant();
  return kelvin*R/1000.;
}

double kelvin2kJpermol(const double kelvin, const Configuration& config) {
  const double R = config.physical_constants().ideal_gas_constant();
  return kelvin*R/1000.;
}

System spce(const argtype& args) {
  System system;
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
    system.add(config);
  }
  system.add(Potential(
    MakeEwald({{"kmax_squared", args_.key("kmax_squared").dflt("38").str()},
               {"alpha",
      str(args_.key("alphaL").dflt("5.6").dble()/
          system.configuration().domain().min_side_length())}})));
  system.add(Potential(MakeModelTwoBodyFactory({MakeLennardJones(),
                                                MakeChargeScreened()})));
  system.add(Potential(MakeChargeScreenedIntra(), MakeVisitModelBond()));
  system.add(Potential(MakeChargeSelf()));
  system.add(Potential(MakeLongRangeCorrections()));
  return system;
}

System rpm(const argtype& args) {
  System system;
  Arguments args_(args);
  double dual_cut = args_.key("dual_cut").dflt("-1").dble();
  const std::string cubic_box_length =
    args_.key("cubic_box_length").dflt("12").str();
  { std::shared_ptr<Domain> domain;
    if (std::abs(dual_cut + 1) < NEAR_ZERO) {
      domain = MakeDomain({{"cubic_box_length", cubic_box_length}});
    } else {
      domain = MakeDomain({{"cubic_box_length", cubic_box_length},
                           {"init_cells", str(dual_cut)}});
    }
    Configuration config(domain,
      {{"particle_type0", install_dir() + "/" + args_.key("particle0").dflt("plugin/ewald/forcefield/data.rpm_plus").str()},
       {"particle_type1", install_dir() + "/" + args_.key("particle1").dflt("plugin/ewald/forcefield/data.rpm_minus").str()}}
    );

    if (args_.key("cutoff").used()) {
      const double cutoff = args_.dble();
      ASSERT(cutoff > dual_cut,
        "cutoff: " << cutoff << " should be > dual_cut: " << dual_cut);
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
    system.add(config);
  }
  system.add(Potential(
    MakeEwald({{"kmax_squared", args_.key("kmax_squared").dflt("38").str()},
               {"alpha",
      str(args_.key("alphaL").dflt("5.6").dble()/
          system.configuration().domain().min_side_length())}})));
  system.add(Potential(MakeModelTwoBodyFactory({MakeHardSphere(),
                                                MakeChargeScreened()})));
  system.add(Potential(MakeChargeSelf()));
//  std::string iref = "-1";
//  std::string num_steps = "1";
  if (std::abs(dual_cut + 1) > NEAR_ZERO) {
    Potential ref(MakeModelTwoBodyFactory({MakeHardSphere(),
                                           MakeChargeScreened()}),
                  MakeVisitModelCell());
    ref.set_model_params(system.configuration());
    ref.set_model_param("cutoff", 0, dual_cut);
    ref.set_model_param("cutoff", 1, dual_cut);
    system.add_to_reference(ref);
//    iref = "0";
//    num_steps = "4";
  }
//  const double temp = args_.key("temperature").dflt("0.047899460618081").dble();
//  system.set(MakeMetropolis({{"beta", str(1./temp)},
//    {"chemical_potential", str(args_.key("beta_mu").dflt("-13.94").dble()/temp)}}));
//  system.add(MakeTrialTranslate({
//    {"weight", "0.25"},
//    {"tunable_param", "0.1"},
//    {"reference_index", iref},
//    {"num_steps", num_steps}}));
  return system;
}

}  // namespace feasst
