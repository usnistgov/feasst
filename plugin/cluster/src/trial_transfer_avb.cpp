#include "utils/include/serialize.h"
#include "monte_carlo/include/perturb_remove.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/perturb_add_avb.h"
#include "cluster/include/compute_add_avb.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/compute_remove_avb.h"
#include "cluster/include/trial_transfer_avb.h"

namespace feasst {

class MapTrialAddAVB {
 public:
  MapTrialAddAVB() {
    auto obj = MakeTrialAddAVB();
    obj->deserialize_map()["TrialAddAVB"] = obj;
  }
};

static MapTrialAddAVB mapper_ = MapTrialAddAVB();

// Note that changes here should also be incorported into TrialGrow
TrialAddAVB::TrialAddAVB(argtype * args) : Trial(args) {
  class_name_ = "TrialAddAVB";
  set_description("TrialAddAVB");
  args->insert({"grand_canonical", "true"});
  auto perturb = std::make_shared<PerturbAddAVB>(args);
  ASSERT(perturb->delay_add(), "ComputeAddAVB assumes delay_add is true");
  add_stage(
    std::make_shared<SelectParticleAVB>(args),
    perturb,
    args);
  set(MakeComputeAddAVB());
}
TrialAddAVB::TrialAddAVB(argtype args) : TrialAddAVB(&args) {
  FEASST_CHECK_ALL_USED(args);
}

TrialAddAVB::TrialAddAVB(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 3649, "mismatch version: " << version);
}

void TrialAddAVB::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(3649, ostr);
}

class MapTrialRemoveAVB {
 public:
  MapTrialRemoveAVB() {
    auto obj = MakeTrialRemoveAVB();
    obj->deserialize_map()["TrialRemoveAVB"] = obj;
  }
};

static MapTrialRemoveAVB mapper_trial_remove_avb_ = MapTrialRemoveAVB();

TrialRemoveAVB::TrialRemoveAVB(argtype * args) : Trial(args) {
  class_name_ = "TrialRemoveAVB";
  set_description("TrialRemoveAVB");
  args->insert({"grand_canonical", "true"});
  add_stage(
    std::make_shared<SelectParticleAVB>(args),
    std::make_shared<PerturbRemove>(),
    args);
  set(MakeComputeRemoveAVB());
}
TrialRemoveAVB::TrialRemoveAVB(argtype args) : TrialRemoveAVB(&args) {
  FEASST_CHECK_ALL_USED(args);
}

TrialRemoveAVB::TrialRemoveAVB(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 3409, "mismatch version: " << version);
}

void TrialRemoveAVB::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(3409, ostr);
}

class MapTrialTransferAVB {
 public:
  MapTrialTransferAVB() {
    auto obj = MakeTrialTransferAVB();
    obj->deserialize_map()["TrialTransferAVB"] = obj;
  }
};

static MapTrialTransferAVB mapper_trial_transfer_avb__ = MapTrialTransferAVB();

TrialTransferAVB::TrialTransferAVB(argtype * args) : TrialFactoryNamed() {
  class_name_ = "TrialTransferAVB";
  argtype orig_args = *args;
  auto trial_add = MakeTrialAddAVB(orig_args);
  trial_add->set_weight(trial_add->weight()/2.);
  add(trial_add);
  auto trial_remove = std::make_shared<TrialRemoveAVB>(args);
  trial_remove->set_weight(trial_remove->weight()/2.);
  add(trial_remove);
}
TrialTransferAVB::TrialTransferAVB(argtype args) : TrialTransferAVB(&args) {
  FEASST_CHECK_ALL_USED(args);
}

}  // namespace feasst
