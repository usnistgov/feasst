#include "utils/include/serialize.h"
#include "utils/include/io.h"
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_select_all.h"
#include "monte_carlo/include/perturb_volume.h"
#include "gibbs/include/compute_gibbs_volume_transfer.h"
#include "gibbs/include/trial_gibbs_volume_transfer.h"

namespace feasst {

class MapTrialGibbsVolumeTransferOneWay {
 public:
  MapTrialGibbsVolumeTransferOneWay() {
    auto obj = MakeTrialGibbsVolumeTransferOneWay({{"to_configuration_index", "1"}});
    obj->deserialize_map()["TrialGibbsVolumeTransferOneWay"] = obj;
  }
};

static MapTrialGibbsVolumeTransferOneWay mapper_trial_remove_avb_ = MapTrialGibbsVolumeTransferOneWay();

TrialGibbsVolumeTransferOneWay::TrialGibbsVolumeTransferOneWay(argtype * args) : Trial(args) {
  class_name_ = "TrialGibbsVolumeTransferOneWay";
  set_description("TrialGibbsVolumeTransferOneWay");
  const int to_configuration_index = integer("to_configuration_index", args);
  const int configuration_index = integer("configuration_index", args, 0);
  const bool uniform_volume = boolean("uniform_volume", args, true);
  ASSERT(uniform_volume,
    "Gibbs volume transfers must be chosen uniformly in V.");
  args->insert({"uniform_volume", str(uniform_volume)});
  argtype args2 = *args;
  args2.insert({"configuration_index", str(to_configuration_index)});
  add_stage(
    std::make_shared<TrialSelectAll>(&args2),
    std::make_shared<PerturbVolume>(&args2),
    &args2);
  feasst_check_all_used(args2);
  args->insert({"configuration_index", str(configuration_index)});
  args->insert({"constrain_volume_change", "true"});
  add_stage(
    std::make_shared<TrialSelectAll>(args),
    std::make_shared<PerturbVolume>(args),
    args);
  set(MakeComputeGibbsVolumeTransfer());
}
TrialGibbsVolumeTransferOneWay::TrialGibbsVolumeTransferOneWay(argtype args) : TrialGibbsVolumeTransferOneWay(&args) {
  feasst_check_all_used(args);
}

TrialGibbsVolumeTransferOneWay::TrialGibbsVolumeTransferOneWay(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 9234, "mismatch version: " << version);
}

void TrialGibbsVolumeTransferOneWay::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(9234, ostr);
}

class MapTrialGibbsVolumeTransfer {
 public:
  MapTrialGibbsVolumeTransfer() {
    auto obj = MakeTrialGibbsVolumeTransfer();
    obj->deserialize_map()["TrialGibbsVolumeTransfer"] = obj;
  }
};

static MapTrialGibbsVolumeTransfer mapper_trial_transfer_avb__ = MapTrialGibbsVolumeTransfer();

TrialGibbsVolumeTransfer::TrialGibbsVolumeTransfer(argtype * args) : TrialFactoryNamed() {
  class_name_ = "TrialGibbsVolumeTransfer";
  const int config0 = integer("configuration_index0", args, 0);
  const int config1 = integer("configuration_index1", args, 1);
  //INFO("args " << str(*args));
  //ASSERT(!used("configuration_index", *args),
  //  "Do not use argument:configuration_index. Use configuration_index0 or 1.");
  argtype args1 = *args;
  args1.insert({"configuration_index", str(config0)});
  args1.insert({"to_configuration_index", str(config1)});
  args->insert({"configuration_index", str(config1)});
  args->insert({"to_configuration_index", str(config0)});
  auto trial1 = MakeTrialGibbsVolumeTransferOneWay(args1);
  trial1->set_weight(trial1->weight()/2.);
  add(trial1);
  auto trial2 = std::make_shared<TrialGibbsVolumeTransferOneWay>(args);
  trial2->set_weight(trial2->weight()/2.);
  add(trial2);
}
TrialGibbsVolumeTransfer::TrialGibbsVolumeTransfer(argtype args) : TrialGibbsVolumeTransfer(&args) {
  feasst_check_all_used(args);
}

}  // namespace feasst
