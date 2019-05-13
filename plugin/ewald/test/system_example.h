
#ifndef FEASST_EWALD_SYSTEM_EXAMPLE_H_
#define FEASST_EWALD_SYSTEM_EXAMPLE_H_

#include "system/include/system.h"
#include "ewald/include/model_charge_self.h"
#include "ewald/include/model_charge_intra.h"
#include "ewald/include/model_charge_screened.h"
#include "system/include/model_lj.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model_intra.h"
#include "configuration/test/configuration_test.h"

namespace feasst {

void add_ewald_with(std::shared_ptr<ModelTwoBody> model,
    Configuration * config,
    PotentialFactory * full) {
  { Potential potential;
    auto factory = std::make_shared<ModelTwoBodyFactory>();
    factory->add_model(model);
    factory->add_model(std::make_shared<ModelChargeScreened>());
    potential.set_model(factory);
    potential.set_visit_model(std::make_shared<VisitModel>());
    full->add_potential(potential);
  }
  { Potential potential;
    potential.set_model(std::make_shared<ModelChargeIntra>());
    auto visitor = std::make_shared<VisitModelIntra>();
    visitor->set_intra_cut(0);
    potential.set_visit_model(visitor);
    full->add_potential(potential);
  }
  { Potential potential;
    potential.set_visit_model(std::make_shared<LongRangeCorrections>());
    full->add_potential(potential);
  }
  { Potential potential;
    potential.set_model(std::make_shared<ModelChargeSelf>());
    potential.set_visit_model(std::make_shared<VisitModel>());
    full->add_potential(potential);
  }
  { Potential potential;
    auto ewald = std::make_shared<Ewald>();
    ewald->set_kmax(27);
    ewald->update_wave_vectors(*config);
    ewald->init_wave_vector_storage(config);
    potential.set_visit_model(ewald);
    full->add_potential(potential);
  }
}

System spce() {
  PotentialFactory full;
  Configuration config = spce_sample();
  config.add_model_param("alpha", 5.6/config.domain().min_side_length());
  add_ewald_with(std::make_shared<ModelLJ>(), &config, &full);
  System sys;
  sys.add(config);
  sys.set_unoptimized(full);
  sys.precompute();
  return sys;
}

}  // namespace feasst

#endif  // FEASST_EWALD_SYSTEM_EXAMPLE_H_
