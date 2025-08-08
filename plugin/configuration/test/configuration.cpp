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
  for (const std::string spce : {
      "../plugin/configuration/test/data/spce.txt.new",
      "../plugin/configuration/test/data/spce.txt.old"}) {
    INFO(spce);
    Configuration config({{"particle_type", "../particle/atom.txt,../particle/lj.txt,../particle/"+spce+",../particle/lj.txt"}});
    //Configuration config(argtype({{"particle_type", "../particle/atom.txt,../particle/lj.txt,../particle/"+spce+",../particle/lj.txt"}}));
    //std::stringstream ss;
    //config.serialize(ss);
    //INFO(ss.str());
  //  TRY(
  //    auto config2 = config;
  //    config2.add_particle_type("../particle/atom.txt");
  //    CATCH_PHRASE("already provided");
  //  );
    EXPECT_EQ(4, config.num_particle_types());
    EXPECT_EQ("../particle/atom.txt", config.type_to_file_name(0));
    EXPECT_EQ("../particle/lj.txt", config.type_to_file_name(1));
    EXPECT_EQ("../particle/" + spce, config.type_to_file_name(2));

    config.add_particle_type("../particle/lj.txt");
    EXPECT_EQ(5, config.num_particle_types());
    EXPECT_EQ("../particle/lj.txt", config.type_to_file_name(4));
    EXPECT_EQ("4", config.particle_type_to_name(4));

    EXPECT_EQ(6, config.num_site_types());
    EXPECT_EQ(0, config.site_type_to_particle_type(0));
    EXPECT_EQ(1, config.site_type_to_particle_type(1));
    EXPECT_EQ(2, config.site_type_to_particle_type(2));
    EXPECT_EQ(2, config.site_type_to_particle_type(3));
    EXPECT_EQ(3, config.site_type_to_particle_type(4));
    EXPECT_EQ(4, config.site_type_to_particle_type(5));
    const std::vector<std::vector<int> > nstppt =
      config.num_site_types_per_particle_type();
    EXPECT_EQ(5, static_cast<int>(nstppt.size()));
    EXPECT_EQ(1, nstppt[0][0]);
    EXPECT_EQ(0, nstppt[1][0]);
    EXPECT_EQ(1, nstppt[1][1]);
    EXPECT_EQ(0, nstppt[2][0]);
    EXPECT_EQ(0, nstppt[2][1]);
    EXPECT_EQ(1, nstppt[2][2]);
    EXPECT_EQ(2, nstppt[2][3]);
    EXPECT_EQ(1, nstppt[3][4]);
    EXPECT_EQ(1, nstppt[4][5]);

    EXPECT_NEAR(config.unique_type(2, 3).position().coord(0), 1., NEAR_ZERO);
    EXPECT_EQ("0", config.site_type_to_name(0));
    EXPECT_EQ("1", config.site_type_to_name(1));
    EXPECT_EQ("1", config.site_index_to_name(1, 0));
    if (spce == "spce.txt" ||
        spce == "../plugin/configuration/test/data/spce.txt.old") {
      EXPECT_EQ("4", config.site_type_to_name(4));
      EXPECT_EQ("5", config.site_type_to_name(5));
    } else {
      EXPECT_EQ("2", config.site_type_to_name(4));
      EXPECT_EQ("3", config.site_type_to_name(5));
    }
    EXPECT_EQ(config.particle_type(0).site(0).name(), "0");
    EXPECT_EQ(config.unique_type(0).site(0).name(), "0");
    EXPECT_EQ(config.particle_type(1).site(0).name(), "1");
    EXPECT_EQ(config.unique_type(1).site(0).name(), "1");
    if (spce == "spce.txt" ||
        spce == "../plugin/configuration/test/data/spce.txt.old") {
      EXPECT_EQ("2", config.site_type_to_name(2));
      EXPECT_EQ("3", config.site_type_to_name(3));
      EXPECT_EQ(config.particle_type(2).site(0).name(), "2");
      EXPECT_EQ(config.particle_type(2).site(1).name(), "3");
      EXPECT_EQ(config.particle_type(2).site(2).name(), "4");
      EXPECT_EQ(config.unique_type(2).site(0).name(), "2");
      EXPECT_EQ(config.unique_type(2).site(1).name(), "3");
      EXPECT_EQ(config.particle_type(2).bond(0).name(), "0");
      EXPECT_EQ(config.particle_type(2).bond(1).name(), "1");
      EXPECT_EQ(config.unique_type(2).bond(0).name(), "0");
      EXPECT_EQ(config.particle_type(2).angle(0).name(), "0");
      EXPECT_EQ(config.unique_type(2).angle(0).name(), "0");
      EXPECT_EQ(config.particle_type(3).site(0).name(), "5");
      EXPECT_EQ(config.particle_type(4).site(0).name(), "6");
    } else if (spce ==  "../plugin/configuration/test/data/spce.txt.new") {
      EXPECT_EQ("O", config.site_type_to_name(2));
      EXPECT_EQ("H", config.site_type_to_name(3));
      EXPECT_EQ("H2", config.site_index_to_name(2, 2));
      EXPECT_EQ(config.particle_type(2).site(0).name(), "O1");
      EXPECT_EQ(config.site_name_to_index("O1"), 0);
      TRY(
        config.site_name_to_index("OX1");
        CATCH_PHRASE("Could not find");
      );
      EXPECT_EQ(config.particle_type(2).site(1).name(), "H1");
      EXPECT_EQ(config.site_name_to_index("H1"), 1);
      EXPECT_EQ(config.particle_type(2).site(2).name(), "H2");
      EXPECT_EQ(config.site_name_to_index("H2"), 2);
      EXPECT_EQ(config.unique_type(2).site(0).name(), "O");
      EXPECT_EQ(config.site_type_name_to_index("O"), 2);
      EXPECT_EQ(config.unique_type(2).site(1).name(), "H");
      EXPECT_EQ(config.site_type_name_to_index("H"), 3);
      EXPECT_EQ(config.unique_type(2).site(1).name(), "H");
      EXPECT_EQ(config.particle_type(2).bond(0).name(), "O1H1");
      EXPECT_EQ(config.particle_type(2).bond(1).name(), "O1H2");
      EXPECT_EQ(config.unique_type(2).bond(0).name(), "OH");
      EXPECT_EQ(config.particle_type(2).angle(0).name(), "H1O1H2");
      EXPECT_EQ(config.unique_type(2).angle(0).name(), "HOH");
      EXPECT_EQ(config.particle_type(3).site(0).name(), "2");
      EXPECT_EQ(config.particle_type(4).site(0).name(), "3");
      EXPECT_EQ(config.site_type_name_to_index("2"), 4);
      EXPECT_EQ(config.site_type_name_to_index("3"), 5);
    }
  }
}

TEST(Configuration, coordinates_and_wrapping) {
  auto config = MakeConfiguration({{"cubic_side_length", "5"},
    {"particle_type", "atom:../particle/atom_new.txt"},
    {"add_num_atom_particles", "2"}});
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
  auto config = MakeConfiguration({{"particle_type", "lj:../particle/lj_new.txt"}});
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
  auto config = MakeConfiguration({{"particle_type", "spce:../particle/spce_new.txt"}});
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
  EXPECT_EQ(0, config->unique_type(0).site(1).property("sigma"));
  EXPECT_EQ(0.4238, config->unique_type(0).site(1).property("charge"));
  EXPECT_EQ(3, config->max_sites_in_any_particle());
}

TEST(Configuration, bonds_spce) {
  Configuration config = spce_sample1();
  EXPECT_EQ(2, config.particle_type(0).num_bonds());
  EXPECT_EQ(1, config.unique_type(0).num_bonds());
  EXPECT_EQ(1., config.unique_type(0).bond(0).property("length"));
  EXPECT_EQ(0.0001, config.unique_type(0).bond(0).property("delta"));
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

TEST(Configuration, order) {
  auto config = MakeConfiguration({{"cubic_side_length", "7"},
                                   {"particle_type", "spce:../particle/spce_new.txt"}});
  config->add_particle_of_type(0);
  TRY(
    config->add_particle_type("../particle/lj.txt");
    CATCH_PHRASE("types cannot be added after particles");
  );
}

TEST(Configuration, group) {
  auto config = MakeConfiguration({{"cubic_side_length", "7"}});
  TRY(
    Configuration config_err(*config);
    config_err.add(MakeGroup({{"site_type", "0"}}));
    CATCH_PHRASE("add groups after particle types");
  );
  config->add_particle_type("../particle/spce_new.txt", "spce");
  TRY(
    Configuration config_err(*config);
    config_err.add(MakeGroup({{"site_type", "Hello"}}));
    CATCH_PHRASE("site_type_name:Hello not found");
  );
  config->add_particle_type("../particle/lj_new.txt", "lj");
  config->add(MakeGroup({{"site_type", "O"}, {"particle_type", "spce"}}), "O");
  config->add(MakeGroup({{"site_type", "LJ"}, {"particle_type", "lj"}}), "LJ");
  config->add(MakeGroup({{"site_type", "H"}, {"particle_type", "spce"}}), "H");
  //config->add(MakeGroup({{"site_type", "2"}, {"particle_type", "1"}}), "none");
  EXPECT_TRUE(config->group_select(0).group().has_property("0"));
  EXPECT_TRUE(config->group_select(1).group().has_property("O"));
  EXPECT_TRUE(config->group_select(2).group().has_property("LJ"));
  EXPECT_TRUE(config->group_select(3).group().has_property("H"));
  //EXPECT_TRUE(config->group_select(3).group().has_property("none"));
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
  EXPECT_EQ(1, sel1.num_particles());
  EXPECT_EQ(1, sel1.num_sites());
  const Select& sel2 = config->group_select(3);
  EXPECT_EQ(100, sel2.num_particles());
  EXPECT_EQ(200, sel2.num_sites());

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
    {"xyz_file", "../plugin/configuration/test/data/lj_sample_config_periodic4.xyz"},
    {"particle_type", "lj:../particle/lj_new.txt"},
    {"group", "first"}, {"first_particle_type", "lj"}});
  EXPECT_EQ(2, config->num_groups());
}

TEST(Configuration, select_particle_by_group) {
  Configuration config = spce_sample1();
  config.add(MakeGroup({{"site_type", "O"}}));
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
  auto config = MakeConfiguration({{"cubic_side_length", "10"},
    {"particle_type", "../particle/spce.txt"}});
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
  auto config = MakeConfiguration({{"particle_type", "lj:../particle/lj.txt,atom:../particle/atom.txt"},
    {"add_num_atom_particles", "2"}});
  EXPECT_EQ(2, config->num_particles());
  EXPECT_EQ(1, config->particle(0).type());
  EXPECT_EQ(1, config->particle(1).type());
}

TEST(Configuration, dihedrals) {
  auto config = MakeConfiguration({{"particle_type", "../particle/n-decane.txt"}, {"add_particles_of_type0", "1"}});
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
  auto config = MakeConfiguration({{"particle_type", "spce:../particle/spce_new.txt"},
    {"sigmaH", "0.1"},
    {"epsilonO", "5.2"},
    {"cutoff", "2.3"},
    {"cutoffO", "12.3"},
    {"chargeH", "3.3"},
    {"chargeH_H", "-7.3"},
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
  auto config = MakeConfiguration({{"particle_type", "lj:../particle/lj_new.txt,atom:../particle/atom_new.txt,spce:../particle/spce_new.txt"}});
  EXPECT_EQ(config->model_params().select("charge").value(0), 0.);
  EXPECT_EQ(config->model_params().select("charge").value(1), 0.);
  EXPECT_EQ(config->model_params().select("charge").value(2), -0.8476);
  EXPECT_EQ(config->model_params().select("charge").value(3), 0.4238);
  config->model_params().check();
}

TEST(Configuration, xyz_euler_file) {
  auto config = MakeConfiguration({
    {"particle_type", "euler:../particle/euler.txt"},
    {"xyz_euler_file", "../plugin/configuration/test/data/euler.xyze"}});
  EXPECT_NEAR(200, config->domain().side_length(0), 1e-6);
  EXPECT_NEAR(35.652456, config->particle(0).site(0).position().coord(0), 1e-6);
  EXPECT_NEAR(-2.7106179, config->particle(0).site(0).euler().phi(), 1e-6);
}

TEST(Configuration, 2dfstprt_with_xyz) {
  auto config = MakeConfiguration({
    {"particle_type", "atom:../particle/atom2d.txt"},
    {"xyz_file", "../plugin/configuration/test/data/2d.xyz"}});
}

TEST(Configuration, param_mixing_file) {
  auto config = MakeConfiguration({
    {"particle_type", "../particle/spce.txt"},
    {"cubic_side_length", "8"},
    {"epsilon_mixing_file", "../plugin/configuration/test/data/eps_mix.txt"}});
  EXPECT_NEAR(1.51, config->model_params().select("epsilon").mixed_value(0, 0), NEAR_ZERO);
  EXPECT_NEAR(6.8, config->model_params().select("epsilon").mixed_value(0, 1), NEAR_ZERO);
}

TEST(Configuration, model_param_file) {
  auto config = MakeConfiguration({
    {"particle_type", "co2:../particle/dimer_mie_CO2.txt,n2:../particle/dimer_mie_N2.txt"},
    {"cubic_side_length", "8"},
    {"model_param_file", "../particle/mie_model_parameters.txt"}
  });
  std::vector<std::string> site_names;
  config->site_type_names(&site_names);
  EXPECT_EQ("CO2", site_names[0]);
  EXPECT_EQ("N", site_names[1]);
  EXPECT_EQ(2, static_cast<int>(site_names.size()));
  EXPECT_NEAR(182.364, config->model_params().select("epsilon").value(0), NEAR_ZERO);
  EXPECT_NEAR(182.364, config->model_params().select("epsilon").mixed_value(0, 0), NEAR_ZERO);
  EXPECT_NEAR(121.5, config->model_params().select("epsilon").mixed_value(0, 1), NEAR_ZERO);
  EXPECT_NEAR(49.590, config->model_params().select("epsilon").mixed_value(1, 1), NEAR_ZERO);
  EXPECT_NEAR(20.27, config->model_params().select("mie_lambda_r").mixed_value(0, 1), NEAR_ZERO);
}

TEST(Configuration, site_name_order) {
  auto config = MakeConfiguration({{"particle_type", "spce:../plugin/configuration/test/data/spce_inverted.txt"},
                                   {"add_num_spce_particles", "1"}});
  EXPECT_EQ(config->unique_type(0).site(0).name(), "H");
  EXPECT_EQ(config->unique_type(0).site(1).name(), "O");
  EXPECT_EQ(config->unique_type(0).site(0).type(), 1);
  EXPECT_EQ(config->unique_type(0).site(1).type(), 0);
  EXPECT_EQ(config->model_params().select("sigma").value(0), 0);
  EXPECT_EQ(config->model_params().select("sigma").value(1), 3.16555789);
  EXPECT_EQ(config->particle_type(0).site(0).name(), "O1");
  EXPECT_EQ(config->particle_type(0).site(1).name(), "H1");
  EXPECT_EQ(config->particle_type(0).site(2).name(), "H2");
  EXPECT_EQ(config->particle_type(0).site(0).type(), 1);
  EXPECT_EQ(config->particle_type(0).site(1).type(), 0);
  EXPECT_EQ(config->particle_type(0).site(2).type(), 0);
  FileXYZ().write("tmp/ttt.xyz", *config);
}

}  // namespace feasst
