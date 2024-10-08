#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "monte_carlo/include/perturb_translate.h"
#include "cluster/include/select_cluster.h"
#include "cluster/include/compute_move_cluster.h"
#include "cluster/include/perturb_rotate_com.h"
#include "cluster/include/trial_translate_cluster.h"

namespace feasst {

FEASST_MAPPER(TrialTranslateCluster,);

TrialTranslateCluster::TrialTranslateCluster(argtype * args) : Trial(args) {
  class_name_ = "TrialTranslateCluster";
  set_description("TrialTranslateCluster");
  add_stage(std::make_shared<SelectCluster>(args),
            std::make_shared<PerturbTranslate>(args));
  set(std::make_shared<ComputeMoveCluster>());
}
TrialTranslateCluster::TrialTranslateCluster(argtype args) : TrialTranslateCluster(&args) {
  feasst_check_all_used(args);
}

TrialTranslateCluster::TrialTranslateCluster(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 4578, "mismatch version: " << version);
}

void TrialTranslateCluster::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(4578, ostr);
}

}  // namespace feasst
