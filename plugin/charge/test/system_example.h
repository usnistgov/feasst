
#ifndef FEASST_EWALD_SYSTEM_EXAMPLE_H_
#define FEASST_EWALD_SYSTEM_EXAMPLE_H_

#include "configuration/test/config_utils.h"
#include "configuration/include/domain.h"
#include "system/include/system.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "system/include/potential.h"
#include "system/include/model_two_body_factory.h"
#include "system/include/visit_model_intra.h"
#include "system/include/visit_model_bond.h"
#include "charge/include/ewald.h"
#include "charge/include/charge_screened.h"
#include "charge/include/charge_screened_intra.h"
#include "charge/include/charge_self.h"

namespace feasst {

inline System chain(const double alpha,
                    const int kmax_squared) {
  System system;
  { auto config = MakeConfiguration({{"cubic_side_length", "20"},
      {"particle_type0", "../particle/chain10_3types.txt"},
      {"particle_type1", "../plugin/charge/particle/rpm_minus.txt"},
      {"particle_type2", "../plugin/charge/particle/rpm_plus.txt"},
      {"add_particles_of_type1", "1"},
      {"add_particles_of_type2", "1"},
    });
    config->update_positions({{0, 0, 0}, {2, 0, 0}});
    double alpha2;
    int kxmax, kymax, kzmax;
    Ewald().tolerance_to_alpha_ks(0.0001, *config, &alpha2, &kxmax, &kymax, &kzmax);
    DEBUG("alpha2 " << alpha2);
//    const double rms = Ewald().fourier_rms(alpha2, 3, config, 0);
//    DEBUG("rms0 " << rms);
    DEBUG("kxmax " << kxmax);
    config->add_model_param("alpha", alpha);
    system.add(config);
  }
  auto ewald= MakeEwald({{"kmax_squared", "27"},
               {"alpha", str(5.6/system.configuration().domain().inscribed_sphere_diameter())}});
  system.add(MakePotential(ewald,
                     {{"prevent_cache", "true"}}));
  system.add(MakePotential(MakeModelTwoBodyFactory(MakeLennardJones(),
                                                   MakeChargeScreened())));
  //system.add(MakePotential(MakeChargeScreenedIntra(),
  //                     MakeVisitModelIntra({{"intra_cut", "0"}})));
  system.add(MakePotential(MakeChargeScreenedIntra(), MakeVisitModelBond()));
  system.add(MakePotential(MakeChargeSelf()));
//  system.add(MakePotential(MakeLongRangeCorrections()));
//  auto ewald = add_ewald_with(MakeLennardJones(), &system, kmax_squared);
  system.add(MakePotential(MakeLennardJones(), MakeVisitModelIntra({{"intra_cut", "1"}})));
//  DEBUG("kxmax " << ewald->kxmax());
//  DEBUG("kymax " << ewald->kymax());
//  DEBUG("kzmax " << ewald->kzmax());
  return system;
}

inline void test_cases(
    /// tuple of constants and expected energy values
    std::vector<std::tuple<std::shared_ptr<PhysicalConstants>, double> > cases,
    std::shared_ptr<Model> model,
    std::shared_ptr<VisitModel> visitor = MakeVisitModel()) {
  for (auto cse : cases) {
    INFO(std::get<0>(cse)->class_name());
    Configuration config = spce_sample1();
    config.add_model_param("alpha", 5.6/config.domain().inscribed_sphere_diameter());
    config.set_physical_constants(std::get<0>(cse));
    Potential potential(model, visitor);
    potential.precompute(&config);
    EXPECT_NEAR(std::get<1>(cse), potential.energy(&config), 1e-10);
    Potential potential2 = test_serialize(potential);
    EXPECT_NEAR(std::get<1>(cse), potential2.energy(&config), 1e-10);
  }
}

}  // namespace feasst

#endif  // FEASST_EWALD_SYSTEM_EXAMPLE_H_
