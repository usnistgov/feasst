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
#include "cluster/include/perturb_add_avb.h"
#include "cluster/include/compute_avb2.h"
#include "cluster/include/compute_avb4.h"
#include "cluster/include/compute_add_avb.h"
#include "cluster/include/compute_remove_avb.h"
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

std::shared_ptr<TrialFactory> MakeTrialGrow(std::vector<argtype> args,
    const argtype& default_args) {
  auto factory = std::make_shared<TrialFactory>(&args[0]);
  // First, determine all trial types from args[0]
  std::vector<std::string> trial_types;
  const int num_args = static_cast<int>(args.size());
  ASSERT(num_args > 0, "TrialGrow requires args.size: " << num_args << " > 0");
  if (used("transfer", args[0])) {
    str("transfer", &args[0]);
    trial_types.push_back("add");
    trial_types.push_back("remove");
  }
  if (used("transfer_avb", args[0])) {
    str("transfer_avb", &args[0]);
    trial_types.push_back("add_avb");
    trial_types.push_back("remove_avb");
  }
  if (used("regrow_avb2", args[0])) {
    str("regrow_avb2", &args[0]);
    trial_types.push_back("regrow_avb2_in");
    trial_types.push_back("regrow_avb2_out");
  }
  if (used("regrow", args[0])) {
    str("regrow", &args[0]);
    trial_types.push_back("regrow");
  }
  if (used("regrow_avb4", args[0])) {
    str("regrow_avb4", &args[0]);
    trial_types.push_back("regrow_avb4");
  }
  if (trial_types.size() == 0) trial_types.push_back("partial_regrow");
  const std::string site = str("site", &args[0], "-1");

  // Second, determine the particle_type and first site from args[0]
  const std::string particle_type = str("particle_type", &args[0]);

  // Finally, add each trial to the factory
  for (const std::string trial_type : trial_types) {
    DEBUG("trial_type: " << trial_type);
    std::shared_ptr<Trial> trial = MakeTrial();
    trial->set_description("TrialGrow" + trial_type);
    std::shared_ptr<TrialCompute> compute = MakeTrialComputeMove();
    for (int iarg = 0; iarg < num_args; ++iarg) {
      DEBUG("iarg: " << iarg);
      argtype iargs = args[iarg];
      std::shared_ptr<TrialSelect> select;
      std::shared_ptr<Perturb> perturb;
      if (iarg == 0 && trial_type != "partial_regrow") {
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
          iargs.insert({"particle_type", particle_type});
          iargs.insert({"site", site});
          iargs.insert({"grand_canonical", "true"});
          select = std::make_shared<SelectParticleAVB>(&iargs);
          auto add_avb = std::make_shared<PerturbAddAVB>(&iargs);
          ASSERT(add_avb->delay_add(), "ComputeAddAVB assumes delay_add");
          perturb = add_avb;
          compute = std::make_shared<ComputeAddAVB>();
        } else if (trial_type == "remove_avb") {
          iargs.insert({"particle_type", particle_type});
          iargs.insert({"site", site});
          iargs.insert({"grand_canonical", "true"});
          select = std::make_shared<SelectParticleAVB>(&iargs);
          perturb = std::make_shared<PerturbRemove>();
          compute = std::make_shared<ComputeRemoveAVB>();
        } else if (trial_type == "regrow_avb2_in" ||
                   trial_type == "regrow_avb2_out") {
          //argtype args_sel, args_mv, args_comp;
          //argtype args_avb2(args[iarg]);
          bool out_to_in;
          if (trial_type == "regrow_avb2_in") {
            out_to_in = true;
            iargs.insert({"out_to_in", "true"});
          } else {
            out_to_in = false;
            iargs.insert({"out_to_in", "false"});
          }
          argtype perturb_args;
          gen_avb2_args_(out_to_in, &iargs, &perturb_args);
          iargs.insert({"particle_type", particle_type});
          iargs.insert({"site", site});
          select = std::make_shared<SelectParticleAVB>(&iargs);
          perturb = std::make_shared<PerturbMoveAVB>(&perturb_args);
          compute = std::make_shared<ComputeAVB2>(&iargs);
        } else if (trial_type == "regrow_avb4") {
          //argtype args_sel, args_mv;
          gen_avb4_args_(&iargs);
          iargs.insert({"particle_type", particle_type});
          iargs.insert({"site", site});
          select = std::make_shared<SelectParticleAVB>(&iargs);
          perturb = std::make_shared<PerturbMoveAVB>(&iargs);
          compute = std::make_shared<ComputeAVB4>();
        } else {
          FATAL("unreocgnized trial_type: " << trial_type);
        }
      } else {
        bool bond = false, angle = false, branch = false;
        if (used("bond", iargs)) {
          str("bond", &iargs);
          bond = true;
        }
        if (used("angle", iargs)) {
          str("angle", &iargs);
          angle = true;
        }
        if (used("branch", iargs)) {
          str("branch", &iargs);
          branch = true;
        }
        ASSERT((bond && !(angle || branch)) ||
               (angle && !(bond || branch)) ||
               (branch && !(bond || angle)), "cannot have two of " <<
          "bond: " << bond << " angle: " << angle << " branch: " << branch);
        if (bond) {
          select = MakeTrialSelectBond({
            {"particle_type", particle_type},
            {"mobile_site", str("mobile_site", &iargs)},
            {"anchor_site", str("anchor_site", &iargs)}});
          perturb = MakePerturbDistance();
        } else if (angle) {
          select = MakeTrialSelectAngle({
            {"particle_type", particle_type},
            {"mobile_site", str("mobile_site", &iargs)},
            {"anchor_site", str("anchor_site", &iargs)},
            {"anchor_site2", str("anchor_site2", &iargs)}});
          perturb = MakePerturbDistanceAngle();
        } else if (branch) {
          select = MakeSelectBranch({
            {"particle_type", particle_type},
            {"mobile_site", str("mobile_site", &iargs)},
            {"mobile_site2", str("mobile_site2", &iargs)},
            {"anchor_site", str("anchor_site", &iargs)},
            {"anchor_site2", str("anchor_site2", &iargs)}});
          perturb = MakePerturbBranch();
        } else {
          FATAL("unrecognized args: " << str(iargs) << ". " <<
                "Requires bond, angle, branch, etc");
        }
      }
      argtype dflt_args = default_args;
      argtype stage_args = {{"num_steps", str("num_steps", &iargs, str("num_steps", &dflt_args, "1"))},
        {"reference_index", str("reference_index", &iargs, str("reference_index", &dflt_args, "-1"))}};
      check_all_used(iargs);
      trial->add_stage(select, perturb, &stage_args);
      check_all_used(stage_args);
    }
    trial->set(compute);
    factory->add(trial);
  }
  return factory;
}

}  // namespace feasst
