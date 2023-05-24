#include "utils/test/utils.h"
#include "configuration/include/configuration.h"
#include "configuration/test/config_utils.h"
#include "configuration/include/file_xyz.h"
#include "configuration/include/domain.h"
#include "utils/include/debug.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "utils/include/debug.h"

namespace feasst {

TEST(Configuration, type_to_file_name) {
  Configuration config({
    {"particle_type0", "../forcefield/atom.fstprt"},
    {"particle_type1", "../forcefield/lj.fstprt"},
    {"particle_type2", "../forcefield/spce.fstprt"},
  });
  TRY(
    auto config2 = config;
    config2.add_particle_type("../forcefield/atom.fstprt");
    CATCH_PHRASE("already provided");
  );
  EXPECT_EQ(3, config.num_particle_types());
  EXPECT_EQ("../forcefield/atom.fstprt", config.type_to_file_name(0));
  EXPECT_EQ("../forcefield/lj.fstprt", config.type_to_file_name(1));
  EXPECT_EQ("../forcefield/spce.fstprt", config.type_to_file_name(2));

  config.add_particle_type("../forcefield/lj.fstprt", "2");
  EXPECT_EQ(4, config.num_particle_types());
  EXPECT_EQ("../forcefield/lj.fstprt2", config.type_to_file_name(3));
}

TEST(Configuration, coordinates_and_wrapping) {
  auto config = MakeConfiguration({{"cubic_box_length", "5"},
    {"particle_type0", "../forcefield/atom.fstprt"},
    {"add_particles_of_type0", "2"}});
  Position pos;
  pos.set_to_origin_3D();
  pos.set_coord(0, -583);
  pos.set_coord(1, 83.34);
  pos.set_coord(2, 0.005783);
  Select select(0, config->select_particle(0));
  //select.particle(0, config);
  config->displace_particles(select, pos);
  config->wrap_particle(0);
  EXPECT_EQ(config->num_particles(), 2);
  EXPECT_EQ(config->dimension(), 3);

  Configuration config2 = test_serialize(*config);

  { const Position& site_position = config2.select_particle(0).site(0).position();
    EXPECT_NEAR(2.,       site_position.coord(0), 10*NEAR_ZERO);
    EXPECT_NEAR(-1.66,    site_position.coord(1), 10*NEAR_ZERO);
    EXPECT_NEAR(0.005783, site_position.coord(2), 10*NEAR_ZERO);
  }

  // test unwrapped
  config2.init_wrap(false);
  config2.displace_particles(select, pos);
  config2.wrap_particle(0);
  Configuration config3 = test_serialize(config2);
  { const Position& site_position = config3.select_particle(0).site(0).position();
    EXPECT_NEAR(2 - 583,       site_position.coord(0), 10*NEAR_ZERO);
    EXPECT_NEAR(-1.66 + 83.34,    site_position.coord(1), 10*NEAR_ZERO);
    EXPECT_NEAR(0.005783 + 0.005783, site_position.coord(2), 10*NEAR_ZERO);
  }
}

TEST(Configuration, particle_types_lj) {
  auto config = MakeConfiguration({{"particle_type0", "../forcefield/lj.fstprt"}});
  config->check();
  EXPECT_EQ(1, config->particle_types().num());
  EXPECT_EQ(1, config->particle_types().num());
  EXPECT_EQ(1, config->particle_types().num_site_types());
  EXPECT_EQ(1, config->particle_types().num_sites());
  EXPECT_EQ(1, config->particle_type(0).num_sites());
  EXPECT_EQ(0, config->particle_type(0).site(0).type());
  EXPECT_EQ(1, config->unique_types().num_site_types());
  EXPECT_EQ(1, config->unique_types().num_sites());
  EXPECT_EQ(0, config->unique_type(0).site(0).type());
  EXPECT_EQ(3., config->unique_type(0).site(0).property("cutoff"));
  EXPECT_EQ(1., config->unique_type(0).site(0).property("epsilon"));
  EXPECT_EQ(1., config->unique_type(0).site(0).property("sigma"));
  EXPECT_EQ(1, config->max_sites_in_any_particle());
}

TEST(Configuration, particle_types_spce) {
  auto config = MakeConfiguration({{"particle_type0", "../forcefield/spce.fstprt"}});
  config->check();
  EXPECT_EQ(2, config->particle_types().num_site_types());
  EXPECT_EQ(3, config->particle_types().num_sites());
  EXPECT_EQ(1, config->particle_types().num());
  EXPECT_EQ(3, config->particle_type(0).num_sites());
  EXPECT_EQ(0, config->particle_type(0).site(0).type());
  EXPECT_EQ(1, config->particle_type(0).site(1).type());
  EXPECT_EQ(1, config->particle_type(0).site(2).type());
  EXPECT_EQ(2, config->unique_types().num_site_types());
  EXPECT_EQ(2, config->unique_types().num_sites());
  EXPECT_EQ(0, config->unique_type(0).site(0).type());
  EXPECT_EQ(1, config->unique_type(0).site(1).type());
  EXPECT_EQ(10., config->unique_type(0).site(0).property("cutoff"));
  EXPECT_EQ(0.650169581, config->unique_type(0).site(0).property("epsilon"));
  EXPECT_EQ(3.16555789, config->unique_type(0).site(0).property("sigma"));
  EXPECT_EQ(-0.8476, config->unique_type(0).site(0).property("charge"));
  EXPECT_EQ(10., config->unique_type(0).site(1).property("cutoff"));
  EXPECT_EQ(0., config->unique_type(0).site(1).property("epsilon"));
  EXPECT_EQ(0., config->unique_type(0).site(1).property("sigma"));
  EXPECT_EQ(0.4238, config->unique_type(0).site(1).property("charge"));
  EXPECT_EQ(3, config->max_sites_in_any_particle());
}

TEST(Configuration, bonds_spce) {
  Configuration config = spce_sample1();
  EXPECT_EQ(2, config.particle_type(0).num_bonds());
  EXPECT_EQ(1, config.unique_type(0).num_bonds());
  EXPECT_EQ(1., config.unique_type(0).bond(0).property("length"));
  EXPECT_EQ(0.000001, config.unique_type(0).bond(0).property("delta"));
  EXPECT_EQ(2, config.particle_type(0).bond(0).num_sites());
  EXPECT_EQ(0, config.particle_type(0).bond(0).site(0));
  EXPECT_EQ(1, config.particle_type(0).bond(0).site(1));
  EXPECT_EQ(2, config.particle_type(0).bond(1).num_sites());
  EXPECT_EQ(0, config.particle_type(0).bond(1).site(0));
  EXPECT_EQ(2, config.particle_type(0).bond(1).site(1));

  // add an actual particle and see if it has more info than necessary
  config.add_particle_of_type(0);
  EXPECT_EQ(0, config.particle(0).num_bonds());
}

TEST(Configuration, group) {
  auto config = MakeConfiguration({{"cubic_box_length", "7"}});
  TRY(
    Configuration config_err(*config);
    config_err.add(MakeGroup({{"site_type", "0"}}));
    CATCH_PHRASE("add groups after particle types");
  );
  config->add_particle_type("../forcefield/spce.fstprt");
  config->add_particle_type("../forcefield/lj.fstprt");
  TRY(
    Configuration config_err(*config);
    config_err.add_particle_of_type(0);
    config_err.add_particle_type("../forcefield/lj.fstprt");
    CATCH_PHRASE("types cannot be added after particles");
  );
  config->add(MakeGroup({{"site_type", "0"}, {"particle_type", "0"}}), "O");
  config->add(MakeGroup({{"site_type", "0"}, {"particle_type", "1"}}), "H");
  config->add(MakeGroup({{"site_type", "2"}, {"particle_type", "1"}}), "none");
  EXPECT_TRUE(config->group_select(0).group().has_property("0"));
  EXPECT_TRUE(config->group_select(1).group().has_property("O"));
  EXPECT_TRUE(config->group_select(2).group().has_property("H"));
  EXPECT_TRUE(config->group_select(3).group().has_property("none"));
  EXPECT_EQ(4, config->group_selects().size());
  for (int part = 0; part < 100; ++part) {
    config->add_particle_of_type(0);
  }
  FileXYZ().load(
    "../plugin/configuration/test/data/spce_sample_config_periodic1.xyz",
    config.get());
  config->add_particle_of_type(1);
  EXPECT_EQ(2, config->num_particle_types());
  EXPECT_EQ(3, config->num_site_types());
  const Select& sel0 = config->group_select(1);
  EXPECT_EQ(100, sel0.num_particles());
  EXPECT_EQ(100, sel0.num_sites());
  for (int index = 0; index < sel0.num_particles(); ++index) {
    EXPECT_EQ(1, sel0.site_indices(index).size());
  }
  const Select& sel1 = config->group_select(2);
  EXPECT_EQ(0, sel1.num_particles());
  EXPECT_EQ(0, sel1.num_sites());
  const Select& sel2 = config->group_select(3);
  EXPECT_EQ(1, sel2.num_particles());
  EXPECT_EQ(1, sel2.num_sites());

//  TrialSelectParticle tselect;
//  tselect.random_particle(config, 1);
//  EXPECT_EQ(1, tselect.mobile().num_particles());
//  EXPECT_EQ(1, tselect.mobile().num_sites());
//  EXPECT_EQ(0, tselect.mobile().site_index(0, 0));
//  EXPECT_GT(100, tselect.mobile().particle_index(0));
//  tselect.random_particle(config, 2);
//  EXPECT_EQ(0, tselect.mobile().num_particles());
//  EXPECT_EQ(0, tselect.mobile().num_sites());
//  tselect.random_particle(config, 3);
//  EXPECT_EQ(1, tselect.mobile().num_particles());
//  EXPECT_EQ(1, tselect.mobile().num_sites());
//  EXPECT_EQ(100, tselect.mobile().particle_index(0));
//
  EXPECT_EQ(301, config->num_sites());
//  config->remove_particles(select);
//  EXPECT_EQ(300, config->num_sites());
}

TEST(Configuration, group_as_arg) {
  auto config = MakeConfiguration({
    {"xyz_file", install_dir() + "/plugin/configuration/test/data/lj_sample_config_periodic4.xyz"},
    {"particle_type0", "../forcefield/lj.fstprt"},
    {"group0", "first"}, {"first_particle_type0", "0"}});
  EXPECT_EQ(2, config->num_groups());
}

TEST(Configuration, select_particle_by_group) {
  Configuration config = spce_sample1();
  config.add(MakeGroup({{"site_type", "0"}}));
  Particle part = config.particle(2, 1);
  EXPECT_EQ(1, part.num_sites());
}

TEST(Configuration, physical_constants) {
  auto config = MakeConfiguration({{"physical_constants", "CODATA2010"}});
  EXPECT_EQ(config->model_params().physical_constants().boltzmann_constant(),
            1.3806488E-23);
  auto config2 = MakeConfiguration();
  EXPECT_EQ(config2->model_params().physical_constants().boltzmann_constant(),
            1.380649E-23);
  TRY(
    auto config3 = MakeConfiguration({{"physical_constants", "bananas"}});
    CATCH_PHRASE("The class name \"bananas\" is not recognized");
  );
}

TEST(Configuration, set_type_type) {
  Configuration config = spce_sample1();
  config.add_particle_of_type(0);
  config.add_particle_of_type(0);
  EXPECT_EQ(1, config.particle_type(0).site(1).type());
  EXPECT_EQ(1, config.particle(0).site(1).type());
  EXPECT_EQ(1, config.particle(1).site(1).type());
  EXPECT_EQ(2, config.particle(1).num_sites_of_type(1));
  config.set_site_type(0, 1, 0);
  EXPECT_EQ(0, config.particle_type(0).site(1).type());
  EXPECT_EQ(0, config.particle(0).site(1).type());
  EXPECT_EQ(0, config.particle(1).site(1).type());
  EXPECT_EQ(1, config.particle(1).num_sites_of_type(1));
}

TEST(Configuration, copy_particles) {
  Configuration config1, config2 = spce_sample1();
  EXPECT_EQ(0, config1.num_particles());
  TRY(
    config1.copy_particles(config2);
    CATCH_PHRASE("does not match");
  );
  config1.copy_particles(config2, true);
  EXPECT_EQ(100, config1.num_particles());
  EXPECT_NEAR(config1.particle(99).site(1).position().coord(1),
    6.158403915960, NEAR_ZERO);
}

TEST(Configuration, change_volume) {
  auto config = MakeConfiguration({{"cubic_box_length", "10"},
    {"particle_type", install_dir() + "/forcefield/spce.fstprt"}});
  config->add_particle_of_type(0);
  config->add_particle_of_type(0);
  Select second;
  second.add_particle(config->particle(1), 1);
  config->displace_particle(second, Position({1, 1, 1}));
  EXPECT_NEAR(config->particle(0).site(0).position().coord(0), 0, NEAR_ZERO);
  EXPECT_NEAR(config->particle(0).site(1).position().coord(0), 1.0, NEAR_ZERO);
  EXPECT_NEAR(config->particle(0).site(2).position().coord(0), -0.333313247568237, NEAR_ZERO);
  EXPECT_NEAR(config->particle(0).site(2).position().coord(1), 0.942816142731718, NEAR_ZERO);
  EXPECT_NEAR(config->particle(1).site(0).position().coord(0), 1, NEAR_ZERO);
  EXPECT_NEAR(config->particle(1).site(1).position().coord(0), 2.0, NEAR_ZERO);
  EXPECT_NEAR(config->particle(1).site(2).position().coord(0), 0.666686752431763, NEAR_ZERO);
  EXPECT_NEAR(config->particle(1).site(2).position().coord(1), 1.942816142731718, NEAR_ZERO);
  EXPECT_EQ(config->domain().volume(), 1000);
  EXPECT_EQ(config->domain().side_length(0), 10);
  EXPECT_EQ(config->domain().side_length(1), 10);
  EXPECT_EQ(config->domain().side_length(2), 10);
  config->change_volume(-10, {{"dimension", "0"}});
  EXPECT_EQ(config->domain().volume(), 990);
  EXPECT_EQ(config->domain().side_length(0), 9.9);
  EXPECT_EQ(config->domain().side_length(1), 10);
  EXPECT_EQ(config->domain().side_length(2), 10);
  EXPECT_NEAR(config->particle(0).site(0).position().coord(0), 0, NEAR_ZERO);
  EXPECT_NEAR(config->particle(0).site(1).position().coord(0), 1.0, NEAR_ZERO);
  EXPECT_NEAR(config->particle(0).site(2).position().coord(0), -0.333313247568237, NEAR_ZERO);
  EXPECT_NEAR(config->particle(0).site(2).position().coord(1), 0.942816142731718, NEAR_ZERO);
  EXPECT_NEAR(config->particle(1).site(0).position().coord(0), 0.99, NEAR_ZERO);
  EXPECT_NEAR(config->particle(1).site(1).position().coord(0), 1.99, NEAR_ZERO);
  EXPECT_NEAR(config->particle(1).site(2).position().coord(0), 0.656686752431763, NEAR_ZERO);
  EXPECT_NEAR(config->particle(1).site(2).position().coord(1), 1.942816142731718, NEAR_ZERO);
  config->change_volume(-10, {{"dimension", "-1"}});
  EXPECT_NEAR(config->domain().volume(), 980, 1e-12);
  EXPECT_NEAR(config->domain().side_length(0), 9.86655379913093, 1e-12);
  EXPECT_NEAR(config->domain().side_length(1), 9.96621595871811, 1e-12);
  EXPECT_NEAR(config->domain().side_length(2), 9.96621595871811, 1e-12);
  EXPECT_NEAR(config->particle(0).site(0).position().coord(0), 0, NEAR_ZERO);
  EXPECT_NEAR(config->particle(0).site(1).position().coord(0), 1.0, NEAR_ZERO);
  EXPECT_NEAR(config->particle(0).site(2).position().coord(0), -0.333313247568237, NEAR_ZERO);
  EXPECT_NEAR(config->particle(0).site(2).position().coord(1), 0.942816142731718, NEAR_ZERO);
  EXPECT_NEAR(config->particle(1).site(0).position().coord(0), 0.986655379913093, NEAR_ZERO);
  EXPECT_NEAR(config->particle(1).site(1).position().coord(0), 1.986655379913093, NEAR_ZERO);
  EXPECT_NEAR(config->particle(1).site(2).position().coord(0), 0.653342132344856, NEAR_ZERO);
  EXPECT_NEAR(config->particle(1).site(2).position().coord(1), 1.93943773860353, 2e-15);
}

TEST(Configuration, add_particles_of_type) {
  auto config = MakeConfiguration({{"particle_type0", install_dir() + "/forcefield/lj.fstprt"},
    {"particle_type1", install_dir() + "/forcefield/atom.fstprt"},
    {"add_particles_of_type1", "2"}});
  EXPECT_EQ(2, config->num_particles());
  EXPECT_EQ(1, config->particle(0).type());
  EXPECT_EQ(1, config->particle(1).type());
}

TEST(Configuration, dihedrals) {
  auto config = MakeConfiguration({{"particle_type", "../forcefield/n-decane.fstprt"}, {"add_particles_of_type0", "1"}});
  EXPECT_EQ(1, config->num_particles());
  EXPECT_EQ(7, config->particle_type(0).num_dihedrals());
  EXPECT_EQ(1, config->unique_type(0).num_dihedrals());
  EXPECT_EQ(0, config->particle(0).num_dihedrals());
  EXPECT_EQ("DihedralTraPPE", config->unique_type(0).dihedral(0).model());
  EXPECT_EQ(8, config->particle_type(0).num_angles());
  EXPECT_EQ(1, config->unique_type(0).num_angles());
  EXPECT_EQ(0, config->particle(0).num_angles());
  EXPECT_EQ("AngleHarmonic", config->unique_type(0).angle(0).model());
  EXPECT_EQ(9, config->particle_type(0).num_bonds());
  EXPECT_EQ(1, config->unique_type(0).num_bonds());
  EXPECT_EQ(0, config->particle(0).num_bonds());
  EXPECT_EQ("RigidBond", config->unique_type(0).bond(0).model());
  EXPECT_EQ(9, config->particle_type(0).num_bonds());
}

TEST(Configuration, set_param) {
  auto config = MakeConfiguration({{"particle_type0", "../forcefield/spce.fstprt"},
    {"sigma1", "0.1"},
    {"epsilon0", "5.2"},
    {"cutoff", "2.3"},
    {"cutoff0", "12.3"},
    {"charge1", "3.3"},
    {"charge1_1", "-7.3"},
    });
  EXPECT_DOUBLE_EQ(config->model_params().select("sigma").value(1), 0.1);
  EXPECT_DOUBLE_EQ(config->model_params().select("epsilon").value(0), 5.2);
  EXPECT_DOUBLE_EQ(config->model_params().select("cutoff").value(0), 12.3);
  EXPECT_DOUBLE_EQ(config->model_params().select("cutoff").value(1), 2.3);
  EXPECT_DOUBLE_EQ(config->model_params().select("charge").value(0), -0.8476);
  EXPECT_DOUBLE_EQ(config->model_params().select("charge").value(1), 3.3);
  EXPECT_DOUBLE_EQ(config->model_params().select("charge").mixed_values()[1][1], -7.3);
}

TEST(Configuration, particle_with_different_params) {
  auto config = MakeConfiguration({
    {"particle_type0", "../forcefield/lj.fstprt"},
    {"particle_type1", "../forcefield/atom.fstprt"},
    {"particle_type2", "../forcefield/spce.fstprt"},
  });
  EXPECT_EQ(config->model_params().select("charge").value(0), 0.);
  EXPECT_EQ(config->model_params().select("charge").value(1), 0.);
  EXPECT_EQ(config->model_params().select("charge").value(2), -0.8476);
  EXPECT_EQ(config->model_params().select("charge").value(3), 0.4238);
  config->model_params().check();
}

TEST(Configuration, xyz_euler_file) {
  auto config = MakeConfiguration({
    {"particle_type0", "../forcefield/euler.fstprt"},
    {"xyz_euler_file", "../plugin/configuration/test/data/euler.xyze"}});
  EXPECT_NEAR(200, config->domain().side_length(0), 1e-6);
  EXPECT_NEAR(35.652456, config->particle(0).site(0).position().coord(0), 1e-6);
  EXPECT_NEAR(-2.7106179, config->particle(0).site(0).euler().phi(), 1e-6);
}

TEST(Configuration, 2dfstprt_with_xyz) {
  auto config = MakeConfiguration({
    //{"particle_type0", "../forcefield/atom2d.fstprt"}});
    {"particle_type0", "../forcefield/atom2d.fstprt"},
    {"xyz_file", "../plugin/configuration/test/data/2d.xyz"}});
//  EXPECT_EQ(2, config->domain().dimension());
//  EXPECT_EQ(2, config->particle_type(0).site(0).position().dimension());
}

}  // namespace feasst
