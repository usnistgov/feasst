#include <gtest/gtest.h>
#include "core/test/particle_test.h"
#include "core/include/particle_factory.h"
#include "core/include/file_lmp.h"

namespace feasst {

TEST(ParticleFactory, size_check) {
  Particle particle = default_particle();
  ParticleFactory particles;
  particles.add(particle);
  Position position;
  position.set_vector({0, 0});
  particle.set_position(position);
  try {
    particles.add(particle);
    particles.check();
    CATCH_PHRASE("size error");
  }
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
  try {
    particles.add("../forcefield/data.atom");
    CATCH_PHRASE("particles by file");
  }
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
}

}  // namespace feasst
