#include "utils/test/utils.h"
#include "math/include/accumulator.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/visit_model_bond.h"
#include "system/include/bond_visitor.h"
#include "system/include/system.h"
#include "system/include/potential.h"
#include "models/include/fene.h"
#include "models/include/lennard_jones_cut_shift.h"

namespace feasst {

System fene_ex() {
  System system;
  system.add(MakeConfiguration({{"cubic_side_length", "12"},
    {"particle_type0", "../plugin/models/particle/chain10.txt"},
    {"add_particles_of_type0", "1"}}));
  // wca
  { auto wca = MakeLennardJonesCutShift();
    ModelParams wca_params = system.configuration().model_params();
    wca->set_wca(0, 0, &wca_params);
    wca->precompute(system.configuration());
    auto potential = MakePotential(wca, MakeVisitModelBond());
    potential->set(wca_params); // use wca_params.
    system.add(potential);
  }
  return system;
}

TEST(FENE, chain10) {
  System system = fene_ex();
  // each bond is at a distance of one
  const double en_per_bond = -0.5*30*1.5*1.5*std::log(1-1/1.5/1.5)+1;
  DEBUG("en_per_bond " << en_per_bond);
  EXPECT_DOUBLE_EQ(9.*en_per_bond, system.energy());
}

double select(Position * xn, Random * ran, FENE * fene, const Bond& bond) {
  int attempt = 0;
  while (attempt < 1e6) {
    ran->position_in_spherical_shell(0, 1.5, xn);
    const double l = xn->distance();
    if (ran->uniform() < std::exp(-1.*fene->energy(l, bond))) {
      return l;
    }
    ++attempt;
  }
  FATAL("max attempts");
}

TEST(FENE, random_LONG) {
  System system = fene_ex();
  Accumulator len1, len2;
  auto ran = MakeRandomMT19937();
  const int num = 1e6;
  auto fene = MakeFENE();
  const Bond& bond = system.configuration().unique_type(0).bond(0);
  for (int i = 0; i < num; ++i) {
    len1.accumulate(fene->random_distance(bond, 1., 3, ran.get()));
  }
  DEBUG(len1.str());

  // unoptimized method
  Position xn;
  xn.set_to_origin(3);
  for (int i = 0; i < num; ++i) {
    len2.accumulate(select(&xn, ran.get(), fene.get(), bond));
  }
  DEBUG(len2.str());
}

}  // namespace feasst
