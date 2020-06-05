#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/domain.h"
#include "system/include/lennard_jones.h"
#include "cluster/include/energy_map_all.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/perturb_add_avb.h"

namespace feasst {

TEST(PerturbMoveAVB, move) {
  System system;
  {
    Configuration config(MakeDomain({{"cubic_box_length", "8"}}),
                         {{"particle_type", "../forcefield/data.lj"}});
    config.add_particle_of_type(0);
    config.add_particle_of_type(0);
    config.update_positions({{0, 0, 0},
                             {1.5, 0, 0}});
    system.add(config);
  }
  const Configuration& config = system.configuration();
  system.add(Potential(MakeLennardJones(),
                       MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
  system.energy();

  auto ran = MakeRandomMT19937();
  // ran = MakeRandomMT19937({{"seed", "default"}});

  const double max_dist = 2;
  auto neigh_crit = MakeNeighborCriteria({{"maximum_distance", str(max_dist)}});

  SelectParticleAVB sel_in(
    neigh_crit,
    {{"grand_canonical", "false"},
     {"particle_type", "0"},
     {"inside", "true"}});
  sel_in.precompute(&system);
  auto sel_in2 = test_serialize(sel_in);

  PerturbMoveAVB mv_out(neigh_crit, {{"inside", "false"}});
  mv_out.precompute(&sel_in2, &system);
  auto mv_out2 = test_serialize(mv_out);

  EXPECT_FALSE(sel_in2.is_ghost());
  EXPECT_TRUE(sel_in2.sel(&system, ran.get()));
  DEBUG(sel_in2.mobile().str());
  DEBUG("pos? " << sel_in2.mobile().particle_positions().size());

  DEBUG("b4 in->out");
  DEBUG(config.particle(0).site(0).position().str());
  DEBUG(config.particle(1).site(0).position().str());
  mv_out2.perturb(&system, &sel_in2, ran.get());
  DEBUG("perturb in->out");
  DEBUG(config.particle(0).site(0).position().str());
  DEBUG(config.particle(1).site(0).position().str());
  mv_out2.revert(&system);
  DEBUG("revert in->out");
  EXPECT_NEAR(config.particle(0).site(0).position().coord(0), 0, NEAR_ZERO);
  EXPECT_NEAR(config.particle(1).site(0).position().coord(0), 1.5, NEAR_ZERO);
  EXPECT_TRUE(sel_in2.sel(&system, ran.get()));
  DEBUG(sel_in2.mobile().str());
  mv_out2.perturb(&system, &sel_in2, ran.get());
  mv_out2.finalize(&system);
  system.energy();
  DEBUG("af in->out");
  DEBUG(config.particle(0).site(0).position().str());
  DEBUG(config.particle(1).site(0).position().str());

  Position rel = config.particle(sel_in2.mobile().particle_index(0)).site(0).position();
  rel.subtract(config.particle(sel_in2.anchor().particle_index(0)).site(0).position());
  EXPECT_FALSE(neigh_crit->is_position_accepted(rel, config.domain()));

  SelectParticleAVB sel_out(
    neigh_crit,
    {{"grand_canonical", "false"},
     {"particle_type", "0"},
     {"inside", "false"}});
  sel_out.precompute(&system);
  auto sel_out2 = test_serialize(sel_out);

  PerturbMoveAVB mv_in(neigh_crit, {{"inside", "true"}});
  mv_in.precompute(&sel_out2, &system);
  auto mv_in2 = test_serialize(mv_in);

  EXPECT_FALSE(sel_out2.is_ghost());
  EXPECT_TRUE(sel_out2.sel(&system, ran.get()));
  DEBUG(sel_out2.mobile().str());

  DEBUG("b4 out->in");
  const double b4pos = config.particle(0).site(0).position().coord(0);
  DEBUG(config.particle(0).site(0).position().str());
  DEBUG(config.particle(1).site(0).position().str());
  mv_in2.perturb(&system, &sel_out2, ran.get());
  DEBUG("b4 revert out->in");
  DEBUG(config.particle(0).site(0).position().str());
  DEBUG(config.particle(1).site(0).position().str());
  mv_in2.revert(&system);
  DEBUG("af revert out->in");
  DEBUG(config.particle(0).site(0).position().str());
  DEBUG(config.particle(1).site(0).position().str());
  EXPECT_NEAR(config.particle(0).site(0).position().coord(0), b4pos, NEAR_ZERO);
  EXPECT_TRUE(sel_out2.sel(&system, ran.get()));
  DEBUG(sel_out2.mobile().str());
  mv_in2.perturb(&system, &sel_out2, ran.get());
  mv_in2.finalize(&system);
  DEBUG("af out->in");
  DEBUG(config.particle(0).site(0).position().str());
  DEBUG(config.particle(1).site(0).position().str());

  rel = config.particle(sel_out2.mobile().particle_index(0)).site(0).position();
  rel.subtract(config.particle(sel_out2.anchor().particle_index(0)).site(0).position());
  EXPECT_TRUE(neigh_crit->is_position_accepted(rel, config.domain()));
}

TEST(PerturbMoveAVB, AVB4) {
  System system;
  {
    Configuration config(MakeDomain({{"cubic_box_length", "8"}}),
                         {{"particle_type", "../forcefield/data.lj"}});
    config.add_particle_of_type(0);
    config.add_particle_of_type(0);
    config.add_particle_of_type(0);
    config.update_positions({{0, 0, 0},
                             {1.5, 0, 0},
                             {-3, 0, 0},
                             });
    system.add(config);
  }
  const Configuration& config = system.configuration();
  system.add(Potential(MakeLennardJones(),
                       MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
  system.energy();

  auto ran = MakeRandomMT19937();
  //ran = MakeRandomMT19937({{"seed", "1580154124"}});
  //ran = MakeRandomMT19937({{"seed", "1588085108"}});

  const double max_dist = 2;
  auto neigh_crit = MakeNeighborCriteria({{"maximum_distance", str(max_dist)}});

  SelectParticleAVB sel_in(
    neigh_crit,
    {{"grand_canonical", "false"},
     {"particle_type", "0"},
     {"inside", "true"},
     {"second_target", "true"}});
  sel_in.precompute(&system);
  auto sel_in2 = test_serialize(sel_in);

  PerturbMoveAVB mv_in(neigh_crit, {{"inside", "true"}});
  mv_in.precompute(&sel_in2, &system);
  auto mv_in2 = test_serialize(mv_in);

  EXPECT_FALSE(sel_in2.is_ghost());
  const bool found = sel_in2.sel(&system, ran.get());
  if (!found) {
    const int mobile = sel_in2.mobile().particle_index(0);
    EXPECT_TRUE(mobile == 0 || mobile == 1);
    return;
  }
  DEBUG(sel_in2.mobile().str());
  DEBUG("pos? " << sel_in2.mobile().particle_positions().size());

  DEBUG("b4");
  DEBUG(config.particle(0).site(0).position().str());
  DEBUG(config.particle(1).site(0).position().str());
  DEBUG(config.particle(2).site(0).position().str());
  mv_in2.perturb(&system, &sel_in2, ran.get());
  mv_in2.finalize(&system);
  system.energy();
  DEBUG("af");
  DEBUG(config.particle(0).site(0).position().str());
  DEBUG(config.particle(1).site(0).position().str());
  DEBUG(config.particle(2).site(0).position().str());
}

}  // namespace feasst
