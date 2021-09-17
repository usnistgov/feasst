#include <string>
#include "utils/test/utils.h"
#include "xtc/include/file_xtc.h"
#include "utils/include/debug.h"

TEST(FileXTC, XTC) {
  feasst::Configuration config;
  feasst::FileXTC file;
  std::string file_name = "../plugin/xtc/test/data/1L2Y.xtc";
  XDRFILE * xdrfile_read_only = xdrfile_open(file_name.c_str(), "r");
  TRY(
    file.load(xdrfile_read_only, &config, "../plugin/xtc/test/data/1L2Y.xtc");
    CATCH_PHRASE("does not match number of atoms in xtc");
  );
  config.add_particle_type("../forcefield/atom.fstprt");
  for (int i = 0; i < 245; ++i) {
    config.add_particle_of_type(0);
  }
  EXPECT_EQ(0, file.load(xdrfile_read_only, &config));
  EXPECT_NEAR(0.8275, config.particle(0).site(0).position().coord(0), 1e-5);
  EXPECT_EQ(1, file.load(xdrfile_read_only, &config));

  XDRFILE * xdrfile_write = xdrfile_open("tmp/1L2Ycpy.xtc", "w");
  file.write(xdrfile_write, config);
}
