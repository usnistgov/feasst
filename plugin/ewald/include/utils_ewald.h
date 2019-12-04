
#ifndef FEASST_EWALD_UTILS_EWALD_H_
#define FEASST_EWALD_UTILS_EWALD_H_

#include "ewald/include/ewald.h"
#include "ewald/include/charge_self.h"
#include "ewald/include/charge_screened_intra.h"
#include "ewald/include/charge_screened.h"
#include "system/include/model_two_body_factory.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model_intra.h"
#include "system/include/system.h"

namespace feasst {

inline std::shared_ptr<Ewald> add_ewald_with(
    std::shared_ptr<ModelTwoBody> model,
    System * system,
    const int kmax_squared = 27) {
  { Potential potential;
    auto factory = MakeModelTwoBodyFactory();
    factory->add_model(model);
    factory->add_model(MakeChargeScreened());
    potential.set_model(factory);
    potential.set_visit_model(MakeVisitModel());
    system->add(potential);
  }
  { Potential potential;
    potential.set_model(MakeChargeScreenedIntra());
    auto visitor = MakeVisitModelIntra();
    visitor->set_intra_cut(0);
    potential.set_visit_model(visitor);
    system->add(potential);
  }
  { Potential potential;
    potential.set_model(MakeChargeSelf());
    potential.set_visit_model(MakeVisitModel());
    system->add(potential);
  }
  auto ewald = MakeEwald();
  { Potential potential;
    ewald->set_kmax_squared(kmax_squared);
    ewald->update_wave_vectors(system->configuration());
    ewald->init_wave_vector_storage(system->get_configuration());
    potential.set_visit_model(ewald);
    system->add(potential);
  }
  return ewald;
}

}  // namespace feasst

#endif  // FEASST_EWALD_UTILS_EWALD_H_
