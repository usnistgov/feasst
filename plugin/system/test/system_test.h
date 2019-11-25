
#ifndef FEASST_CORE_SYSTEM_TEST_H_
#define FEASST_CORE_SYSTEM_TEST_H_

#include <vector>
#include <utility>
#include <string>
#include "system/include/system.h"
#include "system/include/potential.h"
#include "configuration/test/configuration_test.h"
#include "system/include/hard_sphere.h"
#include "system/include/lennard_jones.h"

namespace feasst {

inline Potential default_potential() {
  Potential potential;
  potential.set_model(std::make_shared<LennardJones>());
  potential.set_visit_model(std::make_shared<VisitModel>());
  return potential;
}

inline Potential hs_potential() {
  Potential potential;
  potential.set_model(std::make_shared<HardSphere>());
  potential.set_visit_model(std::make_shared<VisitModel>());
  return potential;
}

inline System default_system() {
  System sys;
  Configuration config = default_configuration();
  sys.add(config);
  sys.add_to_unoptimized(default_potential());
  return sys;
}

}  // namespace feasst

#endif  // FEASST_CORE_SYSTEM_TEST_H_
