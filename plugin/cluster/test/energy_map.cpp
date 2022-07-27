#include "utils/test/utils.h"
#include "utils/include/utils.h"
#include "math/include/constants.h"
#include "math/include/random_mt19937.h"
#include "configuration/test/config_utils.h"
#include "system/include/visit_model.h"
#include "system/include/lennard_jones.h"
#include "cluster/include/energy_map_all.h"
#include "cluster/include/energy_map_all_criteria.h"
#include "cluster/include/energy_map_neighbor.h"
#include "cluster/include/energy_map_neighbor_criteria.h"

namespace feasst {

TEST(EnergyMap, energy_map) {
  const double rcut = 2.;
  auto neighbor_criteria = MakeNeighborCriteria({{"maximum_distance", str(rcut)}});
  //for (std::string mapstr : {"all"}) {
  //for (std::string mapstr : {"all_criteria"}) {
  //for (std::string mapstr : {"neighbor"}) {
  //for (std::string mapstr : {"all", "all_criteria"}) {
  //for (std::string mapstr : {"all", "all_criteria", "neighbor"}) {
  for (std::string mapstr : {"all", "all_criteria", "neighbor", "neighbor_criteria"}) {
    std::shared_ptr<EnergyMap> map;
    if (mapstr == "all") {
      map = MakeEnergyMapAll();
    } else if (mapstr == "all_criteria") {
      map = MakeEnergyMapAllCriteria();
    } else if (mapstr == "neighbor") {
      map = MakeEnergyMapNeighbor();
    } else if (mapstr == "neighbor_criteria") {
      map = MakeEnergyMapNeighborCriteria();
    } else {
      FATAL("unrecognized mapstr");
    }
    Configuration config = lj_sample4();
    config.add(neighbor_criteria);
    LennardJones model;
    model.precompute(config.model_params());
    VisitModel visit(MakeVisitModelInner(map));
    visit.precompute(&config);
    model.compute(&config, &visit);
    visit.finalize(config.selection_of_all(), &config);
    const double en_lj_all = -16.790321304625856;
    EXPECT_NEAR(en_lj_all, visit.energy(), NEAR_ZERO);
    //INFO(visit.inner().energy_map().total_energy());
    if (mapstr == "all" || mapstr == "neighbor") {
      EXPECT_NEAR(en_lj_all,
                  visit.inner().energy_map().total_energy(),
                  1e-13);
    } else if (mapstr == "all_criteria" || mapstr == "neighbor_criteria") {
      EXPECT_NEAR(-15.076312312129398,
                  visit.inner().energy_map().total_energy(),
                  1e-13);
    }
    visit.inner().energy_map().check(config);
    // find neighbors within 3 of first particle manually
    std::vector<int> neighs, neighs_rcut = {1, 5, 8, 11, 14, 25};
    // {1, 4, 5, 8, 9, 11, 13, 14, 15, 21, 25, 27, 28}; // rcut 3
    for (int ipart = 1; ipart < config.num_particles(); ++ipart) {
      Position pos = config.particle(ipart).site(0).position();
      pos.subtract(config.particle(0).site(0).position());
      if (pos.distance() < rcut) {
        neighs.push_back(ipart);
      }
    }
    // INFO("num " << neighs.size());
    // INFO(feasst_str(neighs));
    EXPECT_EQ(neighs, neighs_rcut);

    RandomMT19937 random;
    Select neighs2;
    visit.inner().energy_map().neighbors(
      *neighbor_criteria,
      config,
      0, 0, 0,
      &neighs2);
    //INFO("neighs2 " << neighs2.str());
    const int neighbor = random.const_element(neighs2.particle_indices());
    EXPECT_EQ(neighs_rcut.size(), static_cast<int>(neighs2.num_sites()));
    EXPECT_TRUE(find_in_list(neighbor, neighs2.particle_indices()));

    test_serialize(visit);
  }
}

}  // namespace feasst
