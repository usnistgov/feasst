#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "monte_carlo/include/criteria.h"
#include "steppers/include/criteria_updater.h"

namespace feasst {

CriteriaUpdater::CriteriaUpdater(argtype * args) : ModifyUpdateOnly(args) {}
CriteriaUpdater::CriteriaUpdater(argtype args) : CriteriaUpdater(&args) {
  feasst_check_all_used(args);
}

class MapCriteriaUpdater {
 public:
  MapCriteriaUpdater() {
    CriteriaUpdater().deserialize_map()["CriteriaUpdater"] = MakeCriteriaUpdater();
  }
};

static MapCriteriaUpdater mapper_ = MapCriteriaUpdater();

void CriteriaUpdater::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(743, ostr);
}

CriteriaUpdater::CriteriaUpdater(std::istream& istr) : ModifyUpdateOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 743, "version mismatch:" << version);
}

void CriteriaUpdater::update(Criteria * criteria,
  System * system,
  Random * random,
  TrialFactory * trial_factory) { criteria->update(); }

}  // namespace feasst
