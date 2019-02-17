
#ifndef FEASST_CORE_SYSTEM_TEST_H_
#define FEASST_CORE_SYSTEM_TEST_H_

#include <vector>
#include <utility>
#include <string>
#include "core/include/system.h"
#include "core/include/potential.h"
#include "core/test/configuration_test.h"
#include "core/include/model_hard_sphere.h"
#include "core/include/model_lj.h"

namespace feasst {

inline Potential default_potential() {
  Potential potential;
  potential.set_model(std::make_shared<ModelLJ>());
  potential.set_visit_model(std::make_shared<VisitModel>());
  return potential;
}

inline Potential hs_potential() {
  Potential potential;
  potential.set_model(std::make_shared<ModelHardSphere>());
  potential.set_visit_model(std::make_shared<VisitModel>());
  return potential;
}

inline System default_system() {
  System sys;
  Configuration config = default_configuration();
  sys.add_configuration(config);
  sys.add_to_unoptimized(default_potential());
  return sys;
}

}  // namespace feasst

#endif  // FEASST_CORE_SYSTEM_TEST_H_
