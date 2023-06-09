#include "utils/test/utils.h"
#include "math/include/constants.h"
#include "configuration/test/config_utils.h"
#include "configuration/include/domain.h"
#include "netcdf/include/file_netcdf.h"

namespace feasst {

TEST(FileNETCDF, write) {
  Configuration config = lj_sample4();
  FileNETCDF netcdf({{"file_name", "tmp/print.nc"}, {"float_precision", "false"}});
  //FileNETCDF netcdf("tmp/print.nc", {{"float_precision", "false"}});
  //FileNETCDF netcdf("tmp/print.nc", {{"float_precision", "true"}});
  //FileNETCDF netcdf = FileNETCDF({{"file_name", "tmp/print.nc"}, {"float_precision", "false"}});
  //FileNETCDF netcdf;
  //netcdf = FileNETCDF({{"file_name", "tmp/print.nc"}, {"float_precision", "false"}});
  netcdf.initialize(config);
  netcdf.write(config);
  config.add_particle_of_type(0);
  netcdf.write(config);
  FileNETCDF netcdf2;
}

//TEST(FileNETCDF, load_spce) {
//  Configuration config = spce_sample1();
//  config.check();
//  EXPECT_NEAR(8000, config.domain().volume(), NEAR_ZERO);
//
//  TRY(
//    Configuration conf;
//    FileNETCDF().load("nofile", &conf);
//    CATCH_PHRASE("is empty");
//  );
//}
//
//TEST(FileNETCDF, load_frame) {
//  auto config = MakeConfiguration({{"particle_type0",
//                                    "../forcefield/dimer.fstprt"}});
//  std::ifstream xyz("../plugin/configuration/test/data/dimer4.xyz");
//  FileNETCDF fxyz;
//  fxyz.load_frame(xyz, config.get());
//  EXPECT_EQ(config->particle(0).site(0).position().coord(0), 0);
//  EXPECT_EQ(config->particle(0).site(1).position().coord(0), 0);
//  fxyz.load_frame(xyz, config.get());
//  EXPECT_EQ(config->particle(0).site(0).position().coord(0), 0);
//  EXPECT_EQ(config->particle(0).site(1).position().coord(0), 1);
//  fxyz.load_frame(xyz, config.get());
//  EXPECT_EQ(config->particle(0).site(0).position().coord(0), 0);
//  EXPECT_EQ(config->particle(0).site(1).position().coord(1), 1);
//  fxyz.load_frame(xyz, config.get());
//  EXPECT_EQ(config->particle(0).site(0).position().coord(0), 0);
//  EXPECT_EQ(config->particle(0).site(1).position().coord(2), 1);
//}

}  // namespace feasst
