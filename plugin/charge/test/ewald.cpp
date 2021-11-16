#include <cmath>
#include "utils/test/utils.h"
#include "configuration/test/config_utils.h"
#include "math/include/random_mt19937.h"
#include "system/include/model_empty.h"
#include "system/include/system.h"
#include "charge/include/slab_correction.h"
#include "charge/include/ewald.h"
#include "charge/include/utils.h"
#include "charge/test/charge_utils.h"
#include "charge/test/system_example.h"
#include "monte_carlo/include/perturb_translate.h"
#include "monte_carlo/include/trial_select_particle.h"

namespace feasst {

TEST(Ewald, ewald) {
  Configuration config = spce_sample1();
  auto ewald = MakeEwald({
    {"alpha", str(5.6/config.domain().inscribed_sphere_diameter())},
    {"kmax_squared", "27"}
  });
  ewald->precompute(&config);

  ModelEmpty model;  // any place holder model is fine because its not used
  model.compute(&config, ewald.get());
  ewald->finalize(config.selection_of_all(), &config);
  // ewald.update_eik(config.selection_of_all(), &config);

  //const std::vector<double> eik = config.particle(0).site(0).properties().values();
  const std::vector<double>& eik = ewald->eik()[0][0];
  EXPECT_NEAR(eik[0], 1, NEAR_ZERO);
  EXPECT_NEAR(eik[1], -0.069470287276879206, NEAR_ZERO);
  EXPECT_NEAR(eik[2], -0.99034775837133582, NEAR_ZERO);
  const int ry0_index = 13;
  EXPECT_NEAR(eik[ry0_index + 0], 1, NEAR_ZERO);
  EXPECT_NEAR(eik[ry0_index + 1], -0.87389397051446949, NEAR_ZERO);
  EXPECT_NEAR(eik[ry0_index - 1], eik[ry0_index + 1], NEAR_ZERO);
  const int ix0_index = 33;
  EXPECT_NEAR(eik[ix0_index + 0], 0, NEAR_ZERO);
  EXPECT_NEAR(eik[ix0_index + 1], -0.99758402111584965, NEAR_ZERO);
  EXPECT_NEAR(eik[ix0_index + 2], 0.13860489705948481, NEAR_ZERO);
  const int iz0_index = 59;
  EXPECT_NEAR(eik[iz0_index + 0], 0, NEAR_ZERO);
  EXPECT_NEAR(eik[iz0_index + 1], -0.52837486359383823, NEAR_ZERO);
  EXPECT_NEAR(eik[iz0_index - 1], 0.52837486359383823, NEAR_ZERO);

  EXPECT_NEAR(ewald->struct_fact_real()[0], -1.829963812936731, 5e-15);
  EXPECT_NEAR(ewald->struct_fact_imag()[0], 2.3263016099862206, 5e-15);

  const double en = 52.13245734204947;
  EXPECT_NEAR(ewald->energy(), en, 1e-12);

  EXPECT_NEAR(ewald->net_charge(config), 0., NEAR_ZERO);

  // serialize
  auto ewald2 = test_serialize<Ewald, VisitModel>(*ewald);
  model.compute(&config, ewald2.get());
  ewald2->finalize(Select(), &config);
  EXPECT_NEAR(ewald2->energy(), en, 1e-12);
}

TEST(Ewald, system) {
  const double en_lrc = -6.84874714555147;
  {
    System system = spce({{"physical_constants", "CODATA2010"},
      {"alpha", str(5.6/20)},
      {"kmax_squared", "27"},
      {"xyz_file", "../plugin/configuration/test/data/spce_sample_config_periodic1.xyz"}});
    system.add(MakePotential(MakeSlabCorrection({{"dimension", "0"}})));
    EXPECT_NEAR(-4061.7105206971828, system.energy(), 1e-10);
    EXPECT_NEAR(52.1324574151071, system.potential(0).stored_energy(), 1e-12);
    EXPECT_NEAR(-3819.2497119062596, system.potential(1).stored_energy(), 1e-10);
    EXPECT_NEAR(23363.573774608, system.potential(2).stored_energy(), 1e-10);
    EXPECT_NEAR(-23652.08040365018, system.potential(3).stored_energy(), 1e-12);
    EXPECT_NEAR(en_lrc, system.potential(4).stored_energy(), 1e-13);
    EXPECT_NEAR(0.7621099816966812, system.potential(5).stored_energy(), 1e-12);
  }
  {
    System system = spce({{"physical_constants", "CODATA2018"},
      {"alpha", str(5.6/20)},
      {"kmax_squared", "27"},
      {"xyz_file", "../plugin/configuration/test/data/spce_sample_config_periodic1.xyz"}});
    system.add(MakePotential(MakeSlabCorrection({{"dimension", "0"}})));
    EXPECT_EQ(system.potentials().num(), 6);
    EXPECT_NEAR(-4061.7105138549323, system.energy(), 1e-10);
    EXPECT_NEAR(52.13245734204947, system.potential(0).stored_energy(), 1e-12);
    EXPECT_NEAR(-3819.2497053941993, system.potential(1).stored_energy(), 1e-10);
    EXPECT_NEAR(23363.573741866534, system.potential(2).stored_energy(), 1e-10);
    EXPECT_NEAR(-23652.080370504391, system.potential(3).stored_energy(), 1e-12);
    EXPECT_NEAR(en_lrc, system.potential(4).stored_energy(), 1e-13);
    EXPECT_NEAR(0.76210998062866941, system.potential(5).stored_energy(), 1e-12);
  }
}

TEST(Ewald, revert) {
  System system = spce({{"physical_constants", "CODATA2018"},
    {"alpha", str(5.6/20)},
    {"kmax_squared", "27"},
    {"xyz_file", "../plugin/configuration/test/data/spce_sample_config_periodic1.xyz"}});
  const double en = -4062.4726238355574;
  EXPECT_NEAR(en, system.energy(), 1e-10);
  PerturbTranslate perturb;
  TrialSelectParticle tsel;
  RandomMT19937 random;
  tsel.sel(&system, &random);
  perturb.perturb(&system, &tsel, &random);
  EXPECT_GT(std::abs(en - system.energy()), 1e-10);
  perturb.revert(&system);
  EXPECT_NEAR(en, system.energy(), 1e-10);
}

TEST(Ewald, change_volume) {
  System system = spce({{"alpha", str(5.6/20)}, {"kmax_squared", "38"}});
  try {
    system.change_volume(1);
    CATCH_PHRASE("not implemented");
  }
}

TEST(Ewald, synchronize) {
  System s1 = spce({{"alpha", str(5.6/20)}, {"kmax_squared", "38"}, {"table_size", str(1e3)}});
  s1.get_configuration()->add_particle_of_type(0);
  s1.precompute();
  s1.energy();
  System s2 = test_serialize(s1);
  Select part(0, s1.configuration().particle(0));
  Position disp({0.5, 0.5, 0.5});
  s1.get_configuration()->displace_particle(part, disp);
  s1.energy();

  std::stringstream ss;
  s1.potential(0).visit_model().serialize(ss);
  //INFO(ss.str());
  Ewald ewald1(ss);
//  INFO(s1.potential(0).visit_model().manual_data().dble_3D().size());
//  INFO(ewald1.manual_data().dble_3D().size());
//  std::stringstream ss2;


  EXPECT_NEAR(s1.configuration().particle(0).site(0).position().coord(0), 0.5, NEAR_ZERO);
//  INFO(ewald1.eik().size());
  EXPECT_NEAR(ewald1.eik()[0][0][2], 0.95105651629515364, NEAR_ZERO);
  EXPECT_NEAR(s1.potential(0).visit_model().manual_data().dble_3D()[0][0][2], 0.95105651629515364, NEAR_ZERO);
  EXPECT_NEAR(s2.configuration().particle(0).site(0).position().coord(0), 0., NEAR_ZERO);
  EXPECT_NEAR(s2.potential(0).visit_model().manual_data().dble_3D()[0][0][2], 1, NEAR_ZERO);
  s2.synchronize_(s1, part);
  EXPECT_NEAR(s2.configuration().particle(0).site(0).position().coord(0), 0.5, NEAR_ZERO);
  EXPECT_NEAR(s1.potential(0).visit_model().manual_data().dble_3D()[0][0][2], 0.95105651629515364, NEAR_ZERO);
  EXPECT_NEAR(s2.potential(0).visit_model().manual_data().dble_3D()[0][0][2], 0.95105651629515364, NEAR_ZERO);
}

TEST(Ewald, triclinic) {
  System system = spce({
    //{"xyz_file", "../plugin/charge/test/data/5spce_tilted.xyz"},
    {"xyz_file", "../plugin/charge/test/data/5spce.xyz"},
    {"tolerance", "0.0001"},
  });
  system.energy();
  std::stringstream ss;
  system.potential(0).visit_model().serialize(ss);
  Ewald ewald(ss);
  INFO("alpha " << system.configuration().model_params().property("alpha"));
  INFO("num_vectors " << ewald.num_vectors());
  INFO("num_kx " << ewald.num_kx());
  INFO("num_ky " << ewald.num_ky());
  INFO("num_kz " << ewald.num_kz());
  INFO("kxmax " << ewald.kxmax());
  INFO("kymax " << ewald.kymax());
  INFO("kzmax " << ewald.kzmax());
  INFO("kmax_squared " << ewald.kmax_squared());
  INFO(system.configuration().num_particles());
  INFO(system.stored_energy());
  for (const std::shared_ptr<Potential> pot : system.potentials().potentials()) {
    INFO(pot->visit_model().class_name() << ":" << pot->model().class_name() <<
      " = " << pot->stored_energy());
  }
}

TEST(Ewald, srsw_ref_config) {
  System system;
  system.add(MakeConfiguration({
    {"side_length0", "30.0"},
    {"side_length1", "28.97777478867205"},
    {"side_length2", "29.51512917398008"},
    {"xy", "7.764571353075622"},
    {"xz", "-2.6146722824297473"},
    {"yz", "-4.692615336756641"},
    {"xyz_file", "../plugin/charge/test/data/spce_triclinic_sample_periodic1.xyz"},
    {"particle_type0", install_dir() + "/forcefield/spce.fstprt"},
    {"cutoff", "10"}}));
  system.add(MakePotential(MakeLennardJones()));
  system.add(MakePotential(MakeLongRangeCorrections()));
  system.add(MakePotential(MakeEwald({
    {"alpha", feasst::str(0.2850)},
    {"kxmax", "7"},
    {"kymax", "7"},
    {"kzmax", "7"}})));
  system.add(MakePotential(MakeChargeScreened()));
  system.add(MakePotential(MakeChargeScreenedIntra(), MakeVisitModelBond()));
  system.add(MakePotential(MakeChargeSelf()));
  system.energy();
  EXPECT_NEAR(system.stored_energy_profile()[0], 931.15451, 1e-4);
  EXPECT_NEAR(system.stored_energy_profile()[1], -34.16569, 1e-4);
  EXPECT_NEAR(system.stored_energy_profile()[2], 371.46525, 1e-4);
  EXPECT_NEAR(system.stored_energy_profile()[3], -6046.43627, 1e-4);
  EXPECT_NEAR(system.stored_energy_profile()[4], 95078.89447, 1e-4);
  EXPECT_NEAR(system.stored_energy_profile()[5], -96297.75579, 1e-4);
}

}  // namespace feasst
