#include "utils/test/utils.h"
#include "configuration/include/file_xyz.h"
#include "configuration/test/config_utils.h"
#include "configuration/include/domain.h"
#include "math/include/constants.h"

namespace feasst {

TEST(FileXYZ, load) {
  Configuration config = lj_sample4();
  config.check();
  EXPECT_NEAR(config.domain().volume(), 512., NEAR_ZERO);
  EXPECT_EQ(30, config.num_particles());
  EXPECT_NEAR(config.particle(29).site(0).position().coord(1), 3.786335083587E+00, NEAR_ZERO);
  FileXYZ().write("tmp/print.xyz", config);
}

TEST(FileXYZ, load_spce) {
  Configuration config = spce_sample1();
  config.check();
  EXPECT_NEAR(8000, config.domain().volume(), NEAR_ZERO);

  TRY(
    Configuration conf;
    FileXYZ().load("nofile", &conf);
    CATCH_PHRASE("is empty");
  );
}

TEST(FileXYZ, load_frame) {
  auto config = MakeConfiguration({{"particle_type0",
                                    "../particle/dimer.txt"}});
  std::ifstream xyz("../plugin/configuration/test/data/dimer4.xyz");
  FileXYZ fxyz;
  fxyz.load_frame(xyz, config.get());
  EXPECT_EQ(config->particle(0).site(0).position().coord(0), 0);
  EXPECT_EQ(config->particle(0).site(1).position().coord(0), 0);
  fxyz.load_frame(xyz, config.get());
  EXPECT_EQ(config->particle(0).site(0).position().coord(0), 0);
  EXPECT_EQ(config->particle(0).site(1).position().coord(0), 1);
  fxyz.load_frame(xyz, config.get());
  EXPECT_EQ(config->particle(0).site(0).position().coord(0), 0);
  EXPECT_EQ(config->particle(0).site(1).position().coord(1), 1);
  fxyz.load_frame(xyz, config.get());
  EXPECT_EQ(config->particle(0).site(0).position().coord(0), 0);
  EXPECT_EQ(config->particle(0).site(1).position().coord(2), 1);
}

TEST(FileXYZ, multicomponent) {
  auto config = MakeConfiguration({{"particle_type0", "../particle/lj.txt"},
                                   {"particle_type1", "../particle/spce.txt"},
                                   {"particle_type2", "../particle/lj.txt"},
                                   {"xyz_file", "../plugin/configuration/test/data/multicomponent.xyz"}});
  EXPECT_EQ(config->num_particles(), 3);
  EXPECT_EQ(config->particle(0).type(), 0);
  EXPECT_EQ(config->particle(1).type(), 1);
  EXPECT_EQ(config->particle(2).type(), 2);
  FileXYZ().write("tmp/print.xyz", *config);
}

TEST(FileXYZ, multicomponent_failure) {
  TRY(
    auto config = MakeConfiguration({{"particle_type0", "../particle/lj.txt"},
                                     {"particle_type1", "../particle/spce.txt"},
                                     {"xyz_file", "../plugin/configuration/test/data/multicomponent2.xyz"}});
    CATCH_PHRASE("not given in integers");
  );
}

}  // namespace feasst
