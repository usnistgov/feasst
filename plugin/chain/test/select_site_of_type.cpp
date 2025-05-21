#include "utils/test/utils.h"
#include "monte_carlo/test/monte_carlo_utils.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/configuration.h"
#include "configuration/include/domain.h"
#include "chain/include/select_site_of_type.h"

namespace feasst {

TEST(SelectSiteOfType, serialize) {
  auto config = MakeConfiguration({{"cubic_side_length", "20"},
    {"particle_type", "../particle/chain10_3types.txt"},
    {"add_particles_of_type0", "1"}});
  auto sel = MakeSelectSiteOfType({{"site_type", "0"}});
  Select site;
  auto random = MakeRandomMT19937();
  sel->random_site_in_particle(*config, &site, random.get());
  EXPECT_EQ(site.site_index(0, 0), 0);
  auto sel2 = MakeSelectSiteOfType({{"site_type", "2"}});
  sel2->random_site_in_particle(*config, &site, random.get());
  EXPECT_TRUE(site.site_index(0, 0) == 2 ||
              site.site_index(0, 0) == 9);
  auto sel3 = MakeSelectSiteOfType({{"site_type", "1"}});
  sel3->random_site_in_particle(*config, &site, random.get());
  EXPECT_TRUE(site.site_index(0, 0) == 1 ||
              site.site_index(0, 0) == 3 ||
              site.site_index(0, 0) == 4 ||
              site.site_index(0, 0) == 5 ||
              site.site_index(0, 0) == 6 ||
              site.site_index(0, 0) == 7 ||
              site.site_index(0, 0) == 8);// ||
//              site.site_index(0, 0) == 9);
  SelectSiteOfType sel4 = test_serialize(*sel);
  auto sel5 = MakeSelectSiteOfType({{"site_type", "3"}});
  TRY(
    sel5->random_site_in_particle(*config, &site, random.get());
    CATCH_PHRASE("site_type: 3 is not present");
  );
}

}  // namespace feasst
