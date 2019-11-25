
#ifndef FEASST_FLAT_HISTOGRAM_CONSTRAINT_H_
#define FEASST_FLAT_HISTOGRAM_CONSTRAINT_H_

#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"

namespace feasst {

/**
 */
class Constraint {
 public:
  virtual bool is_allowed(const System* system, const Criteria* criteria) const = 0;
};

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_CONSTRAINT_H_
