#include "utils/test/utils.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/file_particle.h"

namespace feasst {

TEST(ParticleFactory, site_types) {
  ParticleFactory particles;
  particles.add(FileParticle().read("../particle/atom.fstprt"));
  particles.add(FileParticle().read("../particle/spce.fstprt"));
  EXPECT_EQ(2, particles.num_site_types());
  EXPECT_EQ(0, particles.particle(0).site(0).type());
  EXPECT_EQ(0, particles.particle(1).site(0).type());
  EXPECT_EQ(1, particles.particle(1).site(1).type());
  EXPECT_EQ(1, particles.particle(1).site(2).type());
  TRY(
    particles.add("../particle/atom.fstprt");
    CATCH_PHRASE("particles by file");
  );
}

TEST(ParticleFactory, unique_particles) {
  ParticleFactory particles;
  particles.unique_particles();
  particles.add("../particle/atom.fstprt");
  particles.add("../particle/spce.fstprt");
  ParticleFactory part2 = test_serialize(particles);
  EXPECT_EQ(2, part2.num_particle_types());
  EXPECT_EQ(3, part2.num_site_types());
  EXPECT_EQ(1, part2.num_bond_types());
  EXPECT_EQ(1, part2.num_angle_types());
  EXPECT_EQ(0, part2.particle(0).type());
  EXPECT_EQ(0, part2.particle(0).site(0).type());
  EXPECT_EQ(1, part2.particle(1).type());
  EXPECT_EQ(1, part2.particle(1).site(0).type());
  EXPECT_EQ(2, part2.particle(1).site(1).type());
  EXPECT_EQ(2, part2.particle(1).site(2).type());
  EXPECT_EQ(-0.333313247568237000, part2.particle(1).site(2).position().coord(0));
}

}  // namespace feasst
