#include "utils/test/utils.h"
#include "configuration/test/particle_test.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/file_lmp.h"

namespace feasst {

TEST(ParticleFactory, size_check) {
  Particle particle = default_particle();
  ParticleFactory particles;
  particles.add(particle);
  Position position;
  position.set_vector({0, 0});
  particle.set_position(position);
  TRY(
    particles.add(particle);
    particles.check();
    CATCH_PHRASE("size error");
  );
}

TEST(ParticleFactory, site_types) {
  ParticleFactory particles;
  particles.add(FileLMP().read("../forcefield/data.atom"));
  particles.add(FileLMP().read("../forcefield/data.spce"));
  EXPECT_EQ(2, particles.num_site_types());
  EXPECT_EQ(0, particles.particle(0).site(0).type());
  EXPECT_EQ(0, particles.particle(1).site(0).type());
  EXPECT_EQ(1, particles.particle(1).site(1).type());
  EXPECT_EQ(1, particles.particle(1).site(2).type());
  TRY(
    particles.add("../forcefield/data.atom");
    CATCH_PHRASE("particles by file");
  );
}

TEST(ParticleFactory, unique_particles) {
  ParticleFactory particles;
  particles.unique_particles();
  particles.add("../forcefield/data.atom");
  particles.add("../forcefield/data.spce");
  EXPECT_EQ(2, particles.num_particle_types());
  EXPECT_EQ(3, particles.num_site_types());
  EXPECT_EQ(0, particles.particle(0).type());
  EXPECT_EQ(0, particles.particle(0).site(0).type());
  EXPECT_EQ(1, particles.particle(1).type());
  EXPECT_EQ(1, particles.particle(1).site(0).type());
  EXPECT_EQ(2, particles.particle(1).site(1).type());
  EXPECT_EQ(2, particles.particle(1).site(2).type());

  // serialize
  ParticleFactory part2 = test_serialize(particles);
  EXPECT_EQ(3, part2.num_site_types());
  EXPECT_EQ(2, part2.particle(1).site(2).type());
  EXPECT_EQ(-0.333313247568237000, particles.particle(1).site(2).position().coord(0));
}

}  // namespace feasst
