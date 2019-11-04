
#ifndef FEASST_EWALD_UTILS_EWALD_H_
#define FEASST_EWALD_UTILS_EWALD_H_

#include "ewald/include/ewald.h"
#include "ewald/include/model_charge_self.h"
#include "ewald/include/model_charge_intra.h"
#include "ewald/include/model_charge_screened.h"
#include "system/include/model_two_body_factory.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model_intra.h"
#include "system/include/system.h"

namespace feasst {

inline std::shared_ptr<Ewald> add_ewald_with(
    std::shared_ptr<ModelTwoBody> model,
    Configuration * config,
    System * system,
    const int kmax_squared = 27) {
  { Potential potential;
    auto factory = MakeModelTwoBodyFactory();
    factory->add_model(model);
    factory->add_model(MakeModelChargeScreened());
    potential.set_model(factory);
    potential.set_visit_model(MakeVisitModel());
    system->add(potential);
  }
  { Potential potential;
    potential.set_model(MakeModelChargeIntra());
    auto visitor = MakeVisitModelIntra();
    visitor->set_intra_cut(0);
    potential.set_visit_model(visitor);
    system->add(potential);
  }
  { Potential potential;
    potential.set_model(MakeModelChargeSelf());
    potential.set_visit_model(MakeVisitModel());
    system->add(potential);
  }
  auto ewald = MakeEwald();
  { Potential potential;
    ewald->set_kmax_squared(kmax_squared);
    ewald->update_wave_vectors(*config);
    ewald->init_wave_vector_storage(config);
    potential.set_visit_model(ewald);
    system->add(potential);
  }
  return ewald;
}

}  // namespace feasst

#endif  // FEASST_EWALD_UTILS_EWALD_H_
