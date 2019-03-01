
#include <gtest/gtest.h>
#include "core/include/monte_carlo.h"
#include "core/test/monte_carlo_test.h"
#include "flat_histogram/include/criteria_flat_histogram.h"
#include "flat_histogram/include/macrostate_num_particles.h"
#include "flat_histogram/include/bias_wang_landau.h"
#include "core/include/histogram.h"

namespace feasst {

TEST(MonteCarlo, WLMC) {
  MonteCarlo mc = mc_lj();
  mc.add(MakeTrialTransfer({{"weight", "0.25"}}));
  { auto criteria = std::make_shared<CriteriaFlatHistogram>();
    criteria->set_beta(1./1.5);
    criteria->add_activity(exp(-2.775));
    Histogram histogram;
    const int nMol = 5;
    histogram.set_width_center(1., 0.);
    for (int i = 0; i <= nMol; ++i) {
      histogram.add(i);
    }
    auto macrostate = std::make_shared<MacrostateNumParticles>();
    macrostate->set_histogram(histogram);
    criteria->set_macrostate(macrostate);
    auto bias = std::make_shared<BiasWangLandau>();
    bias->resize(histogram);
    criteria->set_bias(bias);
    mc.set(criteria);
  }
  mc.attempt(1e4);
}

}  // namespace feasst
