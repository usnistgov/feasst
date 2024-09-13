#include <cmath>
#include "utils/include/io.h"
#include "utils/include/debug.h"
#include "math/include/constants.h"
#include "configuration/include/domain.h"
#include "configuration/include/file_xyz.h"
#include "system/include/lennard_jones.h"
#include "system/include/hard_sphere.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model_bond.h"
#include "system/include/model_two_body_factory.h"
#include "system/include/visit_model_cell.h"
#include "system/include/visit_model_cutoff_outer.h"
#include "system/include/model_two_body_table.h"
#include "monte_carlo/include/metropolis.h"
#include "charge/test/charge_utils.h"
#include "charge/include/ewald.h"
#include "charge/include/charge_screened.h"
#include "charge/include/charge_screened_intra.h"
#include "charge/include/charge_self.h"
//#include "charge/include/slab_correction.h"
#include "charge/include/trial_add_multiple.h"
#include "charge/include/trial_remove_multiple.h"

namespace feasst {

void add_if_not_used(const std::string& key, argtype * args,
  const std::string& value) {
  if (!used(key, *args)) {
    args->insert({key, value});
  }
}

System spce(argtype args) {
  System system;
  double dual_cut = dble("dual_cut", &args, -1);
  add_if_not_used("cubic_side_length", &args, "20");
  add_if_not_used("physical_constants", &args, "CODATA2018");
  add_if_not_used("particle_type", &args,
    install_dir() + "/particle/spce.fstprt");
  system.add(std::make_shared<Configuration>(&args));
  system.add(MakePotential(std::make_shared<Ewald>(&args)));
  system.add(MakePotential(MakeModelTwoBodyFactory(MakeLennardJones(),
                                                   MakeChargeScreened({{"erfc_table_size", str("erfc_table_size", &args, str(2e4))}})),
                           MakeVisitModelCutoffOuter()//,
                           //{{"table_size", str("table_size", &args, str(1e6))}}
                           ));
  system.add(MakePotential(MakeChargeScreenedIntra(), MakeVisitModelBond()));
  system.add(MakePotential(MakeChargeSelf()));
  if (boolean("lrc", &args, true)) {
    system.add(MakePotential(MakeLongRangeCorrections()));
  }
//  system.add(MakePotential(MakeSlabCorrection({{"dimension", "0"}})));
  if (std::abs(dual_cut + 1) > NEAR_ZERO) {
    std::shared_ptr<Potential> ref;
    ref = MakePotential(MakeModelTwoBodyFactory(MakeLennardJones(),
                                                MakeChargeScreened()),
                        MakeVisitModelCell({{"min_length", str(dual_cut)}}));
    ref->set_model_params(system.configuration());
    ref->set_model_param("cutoff", 0, dual_cut);
    ref->set_model_param("cutoff", 1, dual_cut);
    system.add_to_reference(ref);
  }
  feasst_check_all_used(args);
  return system;
}

System rpm(argtype args) {
  System system;
  double dual_cut = dble("dual_cut", &args, -1);
  add_if_not_used("cubic_side_length", &args, "20");
  add_if_not_used("particle_type0", &args,
    install_dir() + "/plugin/charge/particle/rpm_plus.fstprt");
  add_if_not_used("particle_type1", &args,
    install_dir() + "/plugin/charge/particle/rpm_minus.fstprt");
  {
    Configuration config(&args);
    if (used("cutoff", args)) {
      const double cutoff = dble("cutoff", &args);
      ASSERT(cutoff > dual_cut,
        "cutoff: " << cutoff << " should be > dual_cut: " << dual_cut);
      config.set_model_param("cutoff", 0, cutoff);
      config.set_model_param("cutoff", 1, cutoff);
    }
    const double charge_conversion = config.physical_constants().charge_conversion();
    const ModelParam& charge = config.model_params().select("charge");
    for (int type = 0; type < config.num_particle_types(); ++type) {
      double cval = charge.value(type);
      if (type == 0) {
        cval *= dble("charge_ratio", &args, 1);
      }
      config.set_model_param("charge", type,
        cval/std::sqrt(charge_conversion));
    }
    if (used("delta", args)) {
      const double delta = dble("delta", &args);
      const ModelParam& sigma = config.model_params().select("sigma");
      config.set_model_param("sigma", 0, sigma.value(0) + delta);
      config.set_model_param("sigma", 1, sigma.value(1) - delta);
    }
    system.add(std::make_shared<Configuration>(config));
  }
  system.add(MakePotential(std::make_shared<Ewald>(&args)));
  system.add(MakePotential(MakeModelTwoBodyFactory(MakeHardSphere(),
                                                   MakeChargeScreened())));
  system.add(MakePotential(MakeChargeSelf()));
  if (std::abs(dual_cut + 1) > NEAR_ZERO) {
    auto ref = MakePotential(MakeModelTwoBodyFactory(MakeHardSphere(),
                                                     MakeChargeScreened()),
                        MakeVisitModelCell({{"min_length", str(dual_cut)}}));
    ref->set_model_params(system.configuration());
    ref->set_model_param("cutoff", 0, dual_cut);
    ref->set_model_param("cutoff", 1, dual_cut);
    system.add_to_reference(ref);
  }
  feasst_check_all_used(args);
  return system;
}

}  // namespace feasst
