#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "chain/include/trial_select_site_of_type.h"

namespace feasst {

TEST(TrialSelectSiteOfType, serialize) {
  Configuration config({
    {"particle_type", "../plugin/chain/forcefield/data.chain10titratable"},
    {"cubic_box_length", "20"}});
  config.add_particle_of_type(0);
  auto sel = MakeTrialSelectSiteOfType({{"site_type", "1"}});
  SelectPosition site;
  auto random = MakeRandomMT19937();
  sel->random_site_in_particle(config, &site, random.get());
  EXPECT_EQ(site.site_index(0, 0), 0);
  auto sel2 = MakeTrialSelectSiteOfType({{"site_type", "2"}});
  sel2->random_site_in_particle(config, &site, random.get());
  EXPECT_EQ(site.site_index(0, 0), 2);
  auto sel3 = MakeTrialSelectSiteOfType({{"site_type", "0"}});
  sel3->random_site_in_particle(config, &site, random.get());
  EXPECT_TRUE(site.site_index(0, 0) == 1 ||
              site.site_index(0, 0) == 3 ||
              site.site_index(0, 0) == 4 ||
              site.site_index(0, 0) == 5 ||
              site.site_index(0, 0) == 6 ||
              site.site_index(0, 0) == 7 ||
              site.site_index(0, 0) == 8 ||
              site.site_index(0, 0) == 9);
  TrialSelectSiteOfType sel4 = test_serialize(*sel);
  auto sel5 = MakeTrialSelectSiteOfType({{"site_type", "3"}});
  try {
    sel5->random_site_in_particle(config, &site, random.get());
    CATCH_PHRASE("site_type: 3 is not present");
  }
}

}  // namespace feasst
