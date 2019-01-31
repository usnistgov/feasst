
#ifndef FEASST_CORE_SYSTEM_TEST_H_
#define FEASST_CORE_SYSTEM_TEST_H_

#include <vector>
#include <utility>
#include <string>
#include "core/include/system.h"
#include "core/test/configuration_test.h"
#include "core/include/model_hard_sphere.h"
#include "core/include/model_lj.h"

namespace feasst {

inline Potentials default_potentials() {
  Potential potential;
  potential.set_model(std::make_shared<ModelLJ>());
  potential.set_visit_model(std::make_shared<VisitModel>());
  Potentials potentials;
  potentials.add_potential(potential);
  return potentials;
}

inline Potentials hs_potentials() {
  Potential potential;
  potential.set_model(std::make_shared<ModelHardSphere>());
  potential.set_visit_model(std::make_shared<VisitModel>());
  Potentials potentials;
  potentials.add_potential(potential);
  return potentials;
}

inline System default_system() {
  System sys;
  Configuration config = default_configuration();
  sys.add_configuration(config);
  sys.set_full(default_potentials());
  return sys;
}

}  // namespace feasst

#endif  // FEASST_CORE_SYSTEM_TEST_H_
