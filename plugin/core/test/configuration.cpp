#include <gtest/gtest.h>
#include "core/include/configuration.h"
#include "core/include/file_xyz.h"
#include "core/include/debug.h"

TEST(Configuration, coordinates) {
  feasst::Configuration config;
  config.set_domain(feasst::DomainCuboid().set_cubic(5));
  config.add_particle_type("../forcefield/data.atom");
  config.add_particle(0);
  config.add_particle(0);
  feasst::Position pos;
  pos.set_to_origin_3D();
  pos.set_coord(0, -583);
  pos.set_coord(1, 83.34);
  pos.set_coord(2, 0.005783);
  config.select_particle(1);
  config.displace_selected_particle(pos);
  EXPECT_EQ(config.num_particles(), 2);
  EXPECT_EQ(config.dimension(), 3);
}

TEST(Configuration, particle_types_lj) {
  feasst::Configuration config;
  config.add_particle_type("../forcefield/data.lj");
  EXPECT_EQ(1, config.particle_types().num_site_types());
  EXPECT_EQ(1, config.particle_types().num_sites());
  EXPECT_EQ(0, config.particle_type(0).site(0).type());
  EXPECT_EQ(1, config.unique_types().num_site_types());
  EXPECT_EQ(1, config.unique_types().num_sites());
  EXPECT_EQ(0, config.unique_type(0).site(0).type());
  EXPECT_EQ(3., config.unique_type(0).site(0).property("cutoff"));
  EXPECT_EQ(1., config.unique_type(0).site(0).property("epsilon"));
  EXPECT_EQ(1., config.unique_type(0).site(0).property("sigma"));
}

TEST(Configuration, particle_types_spce) {
  feasst::Configuration config;
  config.add_particle_type("../forcefield/data.spce");
  EXPECT_EQ(2, config.particle_types().num_site_types());
  EXPECT_EQ(3, config.particle_types().num_sites());
  EXPECT_EQ(0, config.particle_type(0).site(0).type());
  EXPECT_EQ(1, config.particle_type(0).site(1).type());
  EXPECT_EQ(1, config.particle_type(0).site(2).type());
  EXPECT_EQ(2, config.unique_types().num_site_types());
  EXPECT_EQ(2, config.unique_types().num_sites());
  EXPECT_EQ(0, config.unique_type(0).site(0).type());
  EXPECT_EQ(1, config.unique_type(0).site(1).type());
  EXPECT_EQ(10., config.unique_type(0).site(0).property("cutoff"));
  EXPECT_EQ(0.650169581, config.unique_type(0).site(0).property("epsilon"));
  EXPECT_EQ(3.16555789, config.unique_type(0).site(0).property("sigma"));
  EXPECT_EQ(-0.8476, config.unique_type(0).site(0).property("charge"));
  EXPECT_EQ(10., config.unique_type(0).site(1).property("cutoff"));
  EXPECT_EQ(0., config.unique_type(0).site(1).property("epsilon"));
  EXPECT_EQ(0., config.unique_type(0).site(1).property("sigma"));
  EXPECT_EQ(0.4238, config.unique_type(0).site(1).property("charge"));

  // bonds
  EXPECT_EQ(2, config.particle_type(0).num_bonds());
  EXPECT_EQ(1, config.unique_type(0).num_bonds());
  EXPECT_EQ(1., config.unique_type(0).bond(0).property("l0"));
  EXPECT_EQ(450., config.unique_type(0).bond(0).property("k"));
  EXPECT_EQ(2, config.particle_type(0).bond(0).num_sites());
  EXPECT_EQ(0, config.particle_type(0).bond(0).site(0));
  EXPECT_EQ(1, config.particle_type(0).bond(0).site(1));
  EXPECT_EQ(2, config.particle_type(0).bond(1).num_sites());
  EXPECT_EQ(0, config.particle_type(0).bond(1).site(0));
  EXPECT_EQ(2, config.particle_type(0).bond(1).site(1));

  // add an actual particle and see if it has more info than necessary
  config.add_particle(0);
  EXPECT_EQ(0, config.particle(0).num_bonds());
}

TEST(Configuration, group) {
  feasst::Configuration config;
  config.set_domain(feasst::DomainCuboid().set_cubic(7));
  try {
    feasst::Configuration config_err(config);
    config_err.add(feasst::Group().add_site_type(0));
    CATCH_PHRASE("add groups after particle types");
  }
  config.add_particle_type("../forcefield/data.spce");
  try {
    feasst::Configuration config_err(config);
    config_err.add_particle(0);
    config_err.add_particle_type("../forcefield/data.lj");
    CATCH_PHRASE("types cannot be added after particles");
  }
  config.add(feasst::Group().add_site_type(0).add_particle_type(0));
  config.add(feasst::Group().add_site_type(0).add_particle_type(1));
  try {
    feasst::Configuration config_err(config);
    config_err.add_particle_type("../forcefield/data.lj");
    CATCH_PHRASE("types cannot be added after groups");
  }
  for (int part = 0; part < 100; ++part) {
    config.add_particle(0);
  }
  feasst::FileXYZ().load("../plugin/core/test/data/spce_sample_config_periodic1.xyz", &config);
  EXPECT_EQ(1, config.num_particle_types());
  EXPECT_EQ(2, config.num_site_types());
  const feasst::Configuration * partial0 = &config.partial(0);
  EXPECT_EQ(1, partial0->num_site_types());
  EXPECT_EQ(-8.384130358330E+00, partial0->particle(0).site(0).position().coord(1));
  EXPECT_EQ(-5.402100413020E+00, partial0->particle(1).site(0).position().coord(1));
  EXPECT_EQ(300, config.num_sites());
  EXPECT_EQ(100, partial0->num_sites());
  EXPECT_EQ(7*7*7, partial0->domain().volume());
  for (std::vector<int> ptf : partial0->partial_to_full_site()) {
    EXPECT_EQ(static_cast<int>(ptf.size()), 1);
    EXPECT_EQ(ptf[0], 0);
  }
  for (std::vector<int> ftps : partial0->full_to_partial_site()) {
    EXPECT_EQ(static_cast<int>(ftps.size()), 3);
    EXPECT_EQ(ftps[0], 0);
    EXPECT_EQ(ftps[1], -1);
    EXPECT_EQ(ftps[2], -1);
  }

  const feasst::Configuration * partial1 = &config.partial(1);
  EXPECT_EQ(0, partial1->num_site_types());
  EXPECT_EQ(0, partial1->num_sites());
  EXPECT_EQ(0, partial1->partial_to_full_site().size());
  EXPECT_EQ(100, partial1->full_to_partial_site().size());
  for (std::vector<int> ftps : partial1->full_to_partial_site()) {
    EXPECT_EQ(static_cast<int>(ftps.size()), 0);
  }
}

TEST(Configuration, cells) {
  feasst::Configuration config;
  config.set_domain(feasst::DomainCuboid().set_cubic(7));
  config.add_particle_type("../forcefield/data.spce");
  config.add(feasst::Group().add_site_type(0));
  config.add_particle(0);
  try {
    config.particle(0).site(0).property("cell");
    CATCH_PHRASE("property not found");
  }
  config.init_cells(1);
  config.init_cells(1.4, 0);
  EXPECT_EQ("cell0", config.domain().cells()[0].label());
  EXPECT_EQ("cell0", config.partial(0).domain().cells()[0].label());
  EXPECT_EQ(config.domain().cells()[0].num_total(), 7*7*7);
  EXPECT_EQ(config.partial(0).domain().cells()[0].num_total(), 5*5*5);
  EXPECT_EQ(7*7*7/2. - 0.5, config.particle(0).site(0).property("cell0"));
  EXPECT_EQ(5*5*5/2. - 0.5, config.partial(0).particle(0).site(0).property("cell0"));
  feasst::Position trajectory({-3.49, -3.49, -3.49});
  config.select_particle(0);
  config.displace_selected_particle(trajectory);
  EXPECT_EQ(0, config.particle(0).site(0).property("cell0"));
  config.select_particle(0);
  config.remove_selected_particle();
  config.check_size();
}

TEST(Configuration, selection) {
  feasst::seed_random_by_date();
  feasst::Configuration config;
  config.set_domain(feasst::DomainCuboid().set_cubic(7));
  config.add_particle_type("../forcefield/data.spce");
  config.add_particle_type("../forcefield/data.lj");
  config.add_particle(0);
  config.add_particle(0);
  config.add_particle(0);
  config.select_random_particle();
  EXPECT_EQ(3, config.num_particles());
  EXPECT_EQ(1, config.selection().num_particles());
  feasst::Selection select = config.selection();
  config.remove_selected_particle();
  EXPECT_EQ(2, config.num_particles());
  EXPECT_EQ(0, config.selection().num_particles());
  try {
    config.set_selection(select);
    CATCH_PHRASE("an expired selection");
  }
  config.select_random_particle_of_type(1);
  EXPECT_TRUE(config.selection().empty());
  config.select_random_particle_of_type(0);
  EXPECT_EQ(1, config.selection().num_particles());
  try {
    config.select_random_particle_of_type(10);
    CATCH_PHRASE("doesn't exist");
  }

  // test select_random_particle in config by histogram of visisted particle from config with more than one particle type (group?) or site.
  config.add_particle(1);
  config.add_particle(1);
  config.add_particle(1);
  EXPECT_EQ(5, config.num_particles());
  int sum = 0;
  int num = 20;
  for (int i = 0; i < num; ++i) {
    config.select_random_particle(feasst::Group().add_particle_type(0));
    int ipart = config.selection().selection()[0].first;
    sum += ipart;
    EXPECT_GE(1, ipart);
  }
  EXPECT_GT(num, sum);
}

TEST(Configuration, displace_selection) {
  feasst::Configuration config;
  config.set_domain(feasst::DomainCuboid().set_cubic(10));
  config.add_particle_type("../forcefield/data.chain10");
  config.add_particle(0);
  config.select_site(0, 0);
  feasst::Position disp;
  disp.set_to_origin_3D();
  disp.set_coord(0, 1.234);
  EXPECT_EQ(0., config.particle(0).site(0).position().coord(0));
  config.displace_selection(disp);
  EXPECT_EQ(1.234, config.particle(0).site(0).position().coord(0));
  config.select_site(0, 9);
  config.displace_selection(disp);
  EXPECT_EQ(1.234+9, config.particle(0).site(9).position().coord(0));
  EXPECT_EQ(8, config.particle(0).site(8).position().coord(0));
}

