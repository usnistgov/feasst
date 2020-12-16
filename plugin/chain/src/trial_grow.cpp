#include "utils/include/debug.h"
#include "utils/include/io.h"
#include "monte_carlo/include/trial_move.h"
#include "monte_carlo/include/trials.h"
#include "monte_carlo/include/trial_compute_add.h"
#include "monte_carlo/include/trial_compute_remove.h"
#include "monte_carlo/include/trial_compute_move.h"
#include "monte_carlo/include/perturb_add.h"
#include "monte_carlo/include/perturb_remove.h"
#include "monte_carlo/include/perturb_anywhere.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/trial_select_bond.h"
#include "monte_carlo/include/trial_select_angle.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/perturb_move_avb.h"
#include "cluster/include/compute_avb2.h"
#include "cluster/include/compute_avb4.h"
#include "cluster/include/trial_avb2.h"
#include "cluster/include/trial_avb4.h"
#include "chain/include/select_end_segment.h"
#include "chain/include/perturb_pivot.h"
#include "chain/include/select_segment.h"
#include "chain/include/perturb_crankshaft.h"
#include "chain/include/select_reptate.h"
#include "chain/include/perturb_reptate.h"
#include "chain/include/select_site_of_type.h"
#include "chain/include/perturb_site_type.h"
#include "chain/include/trial_grow.h"
#include "chain/include/select_branch.h"
#include "chain/include/perturb_branch.h"

namespace feasst {

std::shared_ptr<TrialFactory> MakeTrialGrow(
    const std::vector<argtype>& args,
    const argtype& default_args) {
  auto factory = std::make_shared<TrialFactory>(args[0]);
  // First, determine all trial types from args[0]
  std::vector<std::string> trial_types;
  const int num_args = static_cast<int>(args.size());
  ASSERT(num_args > 0, "TrialGrow requires args.size: " << num_args << " >0");
  Arguments args0(args[0]);
  args0.dont_check();
  if (args0.key("transfer").used()) {
    trial_types.push_back("add");
    trial_types.push_back("remove");
  }
  if (args0.key("transfer_avb").used()) {
    trial_types.push_back("add_avb");
    trial_types.push_back("remove_avb");
  }
  if (args0.key("regrow_avb2").used()) {
    trial_types.push_back("regrow_avb2_in");
    trial_types.push_back("regrow_avb2_out");
  }
  if (args0.key("regrow").used()) trial_types.push_back("regrow");
  if (args0.key("regrow_avb4").used()) trial_types.push_back("regrow_avb4");
  if (trial_types.size() == 0) trial_types.push_back("partial_regrow");
  Arguments default_args_(default_args);
  default_args_.dont_check();

  // Second, determine the particle_type and first site from args[0]
  const std::string particle_type = args0.key("particle_type").str();

  // Finally, add each trial to the factory
  for (const std::string trial_type : trial_types) {
    DEBUG("trial_type: " << trial_type);
    std::shared_ptr<Trial> trial = MakeTrial();
    trial->set_description("TrialGrow" + trial_type);
    std::shared_ptr<TrialCompute> compute = MakeTrialComputeMove();
    for (int iarg = 0; iarg < num_args; ++iarg) {
      DEBUG("iarg: " << iarg);
      Arguments args_(args[iarg]);
      args_.dont_check();
      argtype sel_args;
      std::shared_ptr<TrialSelect> select;
      std::shared_ptr<Perturb> perturb;
      if (iarg == 0 && trial_type != "partial_regrow") {
        const std::string site = args0.key("site").str();
        argtype sel_args = {{"particle_type", particle_type}, {"site", site}};
        select = MakeTrialSelectParticle(sel_args);
        if (trial_type == "add") {
          perturb = MakePerturbAdd();
          compute = MakeTrialComputeAdd();
        } else if (trial_type == "remove") {
          perturb = MakePerturbRemove();
          compute = MakeTrialComputeRemove();
        } else if (trial_type == "regrow") {
          perturb = MakePerturbAnywhere();
        } else if (trial_type == "add_avb") {
          FATAL("add_avb not impl");
        } else if (trial_type == "remove_avb") {
          FATAL("remove avb not impl");
        } else if (trial_type == "regrow_avb2_in" ||
                   trial_type == "regrow_avb2_out") {
          argtype args_sel, args_mv, args_comp;
          argtype args_avb2(args[iarg]);
          if (trial_type == "regrow_avb2_in") {
            args_avb2.insert({"out_to_in", "true"});
          } else {
            args_avb2.insert({"out_to_in", "false"});
          }
          gen_avb2_args_(args_avb2, &args_sel, &args_mv, &args_comp);
          //INFO(str(args_sel));
          select = MakeSelectParticleAVB(args_sel);
          perturb = MakePerturbMoveAVB(args_mv);
          compute = MakeComputeAVB2(args_comp);
        } else if (trial_type == "regrow_avb4") {
          argtype args_sel, args_mv;
          gen_avb4_args_(args[iarg], &args_sel, &args_mv);
          select = MakeSelectParticleAVB(args_sel);
          perturb = MakePerturbMoveAVB(args_mv);
          compute = MakeComputeAVB4();
        } else {
          FATAL("unreocgnized trial_type: " << trial_type);
        }
      } else {
        bool bond = false, angle = false, branch = false;
        if (args_.key("bond").used()) bond = true;
        if (args_.key("angle").used()) angle = true;
        if (args_.key("branch").used()) branch = true;
        ASSERT((bond && !(angle || branch)) ||
               (angle && !(bond || branch)) ||
               (branch && !(bond || angle)), "cannot have two of " <<
          "bond: " << bond << " angle: " << angle << " branch: " << branch);
        if (bond) {
          select = MakeTrialSelectBond({
            {"particle_type", particle_type},
            {"mobile_site", args_.key("mobile_site").str()},
            {"anchor_site", args_.key("anchor_site").str()}});
          perturb = MakePerturbDistance();
        } else if (angle) {
          select = MakeTrialSelectAngle({
            {"particle_type", particle_type},
            {"mobile_site", args_.key("mobile_site").str()},
            {"anchor_site", args_.key("anchor_site").str()},
            {"anchor_site2", args_.key("anchor_site2").str()}});
          perturb = MakePerturbDistanceAngle();
        } else if (branch) {
          select = MakeSelectBranch({
            {"particle_type", particle_type},
            {"mobile_site", args_.key("mobile_site").str()},
            {"mobile_site2", args_.key("mobile_site2").str()},
            {"anchor_site", args_.key("anchor_site").str()},
            {"anchor_site2", args_.key("anchor_site2").str()}});
          perturb = MakePerturbBranch();
        } else {
          FATAL("unrecognized args: " << args_.str() << ". " <<
                "Requires bond, angle, branch, etc");
        }
      }
      trial->add_stage(select, perturb, {
        {"num_steps", args_.key("num_steps").dflt(default_args_.key("num_steps").dflt("1").str()).str()},
        {"reference_index", args_.key("reference_index").dflt(default_args_.key("reference_index").dflt("-1").str()).str()}});
    }
    trial->set(compute);
    factory->add(trial);
  }
  return factory;
}

}  // namespace feasst
