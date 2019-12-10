
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
  system->add(Potential(MakeModelTwoBodyFactory({model, MakeChargeScreened()})));
  system->add(Potential(MakeChargeScreenedIntra(),
                        MakeVisitModelIntra({{"cutoff", "0"}})));
  system->add(Potential(MakeChargeSelf()));
  auto ewald = MakeEwald();
  { ewald->set_kmax_squared(kmax_squared);
    ewald->update_wave_vectors(system->configuration());
    ewald->init_wave_vector_storage(system->get_configuration());
    system->add(Potential(ewald));
  }
  return ewald;
}

}  // namespace feasst

#endif  // FEASST_EWALD_UTILS_EWALD_H_
