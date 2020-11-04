#include <algorithm>  // is_sorted
#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_add.h"
#include "ewald/include/trial_add_multiple.h"
#include "ewald/include/compute_add_multiple.h"

namespace feasst {

std::shared_ptr<Trial> MakeTrialAddMultiple(
    const argtype &args) {
  auto trial = MakeTrial(args);
  trial->set_description("TrialAddMultiple");
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
    trial->add_stage(
      MakeTrialSelectParticle(arg),
      MakePerturbAdd(arg),
      arg);
  }
  trial->set(std::make_shared<ComputeAddMultiple>());
  return trial;
}

std::vector<int> ptypes(Arguments * args) {
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
  DEBUG("ptypes " << feasst_str(ptypes));
  ASSERT(std::is_sorted(ptypes.begin(), ptypes.end()),
    "ptypes not sorted: " << feasst_str(ptypes));
  return ptypes;
}

}  // namespace feasst
