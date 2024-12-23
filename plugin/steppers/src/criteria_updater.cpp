#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/criteria_updater.h"

namespace feasst {

CriteriaUpdater::CriteriaUpdater(argtype * args) : ModifyUpdateOnly(args) {}
CriteriaUpdater::CriteriaUpdater(argtype args) : CriteriaUpdater(&args) {
  feasst_check_all_used(args);
}

FEASST_MAPPER(CriteriaUpdater,);

void CriteriaUpdater::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(743, ostr);
}

CriteriaUpdater::CriteriaUpdater(std::istream& istr) : ModifyUpdateOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 743, "version mismatch:" << version);
}

void CriteriaUpdater::update(MonteCarlo * mc) {
  mc->get_criteria()->update();
}

}  // namespace feasst
