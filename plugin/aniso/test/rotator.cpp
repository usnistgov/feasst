#include "utils/test/utils.h"
#include "system/include/hard_sphere.h"
#include "system/include/lennard_jones.h"
#include "system/include/potential.h"
#include "aniso/include/rotator.h"

namespace feasst {

double contact_dist_rot(std::shared_ptr<ModelTwoBody> model, const double hard_limit_u) {
  System sys;
  sys.add(std::make_shared<Configuration>(argtype({{"cubic_side_length", "100"},
    {"particle_type", "atom:../particle/atom.txt"}, {"add_num_atom_particles", "2"}})));
  sys.add(MakePotential(model));
  Rotator rot({{"contact_tolerance", "1e-8"}, {"hard_limit_u", str(hard_limit_u)}});
  rot.init(&sys, "", "");
  rot.gen_orientations(1, sys.configuration());
  EXPECT_EQ(rot.num_orientations(), 4);
  rot.unique_[0] = -1;
  return rot.contact_distance(0, &sys);
}

TEST(Rotator, contact_distance_lj) {
  const double rh = 0.75;
  EXPECT_NEAR(contact_dist_rot(MakeLennardJones(), 4*(std::pow(rh, -12) - std::pow(rh, -6))), rh, 1e-7);
}

TEST(Rotator, contact_distance_hs) {
  EXPECT_NEAR(contact_dist_rot(MakeHardSphere(), 1e10), 1, 1e-7);
}

TEST(Rotator, unique) {
  System sys;
  sys.add(std::make_shared<Configuration>(argtype({{"cubic_side_length", "100"},
    {"particle_type", "trimer:../particle/trimer.txt"}, {"add_num_trimer_particles", "2"}})));
  sys.add(MakePotential(MakeLennardJones()));
  Rotator rot({{"contact_tolerance", "1e-8"}, {"hard_limit_u", "100"}});
  rot.init(&sys, "", "");
  rot.gen_unique_orientations(1, &sys);
  EXPECT_EQ(rot.num_orientations(), 72);
  EXPECT_EQ(rot.fraction_unique(), 8./72.);
  //rot.unique_[0] = -1;
  //rot.contact_distance(0, &sys);
}

//TEST(Rotator, bounds) {
//  System sys;
//  sys.add(std::make_shared<Configuration>(argtype({{"cubic_side_length", "20"},
//    {"particle_type", "trimer:../particle/trimer.txt"}, {"add_num_trimer_particles", "2"}})));
//  sys.add(MakePotential(MakeLennardJones()));
//  Rotator rot({{"contact_tolerance", "1e-8"}, {"hard_limit_u", "100"}});
//  rot.init(&sys, "", "");
//  //std::vector<std::vector<double> > bounds = {};
//  std::vector<std::vector<double> > bounds = {{0, PI/2.}, {0, PI/2.}, {0, PI/2.}, {0, PI/2.}, {0, PI/2.}};
//  rot.gen_unique_orientations(2, &sys, bounds);
//  //INFO("uni " << rot.num_unique() << " frac " << rot.fraction_unique());
//}

}  // namespace feasst
