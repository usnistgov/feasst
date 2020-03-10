#include <algorithm>  // is_sorted
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_add.h"
#include "ewald/include/trial_add_multiple.h"
#include "ewald/include/compute_add_multiple.h"

namespace feasst {

TrialAddMultiple::TrialAddMultiple(const argtype& args) : Trial(args) {
  Arguments args_(args);
  args_.dont_check();
  // Create one stage per particle
  const std::vector<int> pt = ptypes(&args_);
  std::vector<argtype> new_args;
  for (int p : pt) {
    argtype nag = args_.args();
    nag.insert({"particle_type", str(p)});
    nag.insert({"exclude_perturbed", "true"});
    new_args.push_back(nag);
  }
  for (const argtype& arg : new_args) {
    add_stage(
      MakeTrialSelectParticle(arg),
      MakePerturbAdd(arg),
      arg);
  }
  set(std::make_shared<ComputeAddMultiple>());
  class_name_ = "TrialAddMultiple";
}

std::vector<int> TrialAddMultiple::ptypes(Arguments * args) const {
  std::vector<int> ptypes;
  int count = 0;
  std::string start("particle_type");
  std::stringstream ss;
  ss << start << count;
  while (args->key(ss.str()).used()) {
    DEBUG("ss " << ss.str());
    ptypes.push_back(args->remove().integer());
    ASSERT(count < 1e8, "count: " << count << " is too high");
    ++count;
    ss.str("");
    ss << start << count;
  }
  ASSERT(std::is_sorted(ptypes.begin(), ptypes.end()),
    "ptypes not sorted: " << feasst_str(ptypes));
  DEBUG("ptypes " << feasst_str(ptypes));
  return ptypes;
}

class MapTrialAddMultiple {
 public:
  MapTrialAddMultiple() {
    auto obj = MakeTrialAddMultiple({{"particle_type0", "0"},
                                 {"particle_type1", "1"}});
    obj->deserialize_map()["TrialAddMultiple"] = obj;
  }
};

static MapTrialAddMultiple mapper_ = MapTrialAddMultiple();

std::shared_ptr<Trial> TrialAddMultiple::create(std::istream& istr) const {
  return std::make_shared<TrialAddMultiple>(istr);
}

TrialAddMultiple::TrialAddMultiple(std::istream& istr) : Trial(istr) {
  // ASSERT(class_name_ == "TrialAddMultiple", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(736 == version, "mismatch version: " << version);
}

void TrialAddMultiple::serialize_trial_add_multiple_(std::ostream& ostr) const {
  serialize_trial_(ostr);
  feasst_serialize_version(736, ostr);
}

void TrialAddMultiple::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_add_multiple_(ostr);
}

}  // namespace feasst
