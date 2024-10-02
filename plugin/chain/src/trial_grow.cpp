#include "utils/include/debug.h"
#include "utils/include/arguments_extra.h"
#include "utils/include/serialize_extra.h" // deep_copy
#include "utils/include/io.h"
#include "utils/include/file.h"
#include "monte_carlo/include/trial_move.h"
#include "monte_carlo/include/trial_compute_add.h"
#include "monte_carlo/include/trial_compute_remove.h"
#include "monte_carlo/include/trial_compute_move.h"
#include "monte_carlo/include/perturb_add.h"
#include "monte_carlo/include/perturb_remove.h"
#include "monte_carlo/include/perturb_anywhere.h"
#include "monte_carlo/include/perturb_distance.h"
#include "monte_carlo/include/perturb_distance_angle.h"
#include "monte_carlo/include/perturb_dihedral.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/trial_select_bond.h"
#include "monte_carlo/include/trial_select_angle.h"
#include "monte_carlo/include/trial_select_dihedral.h"
#include "monte_carlo/include/trial_compute_translate.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/perturb_move_avb.h"
#include "cluster/include/perturb_add_avb.h"
#include "cluster/include/compute_avb2.h"
#include "cluster/include/compute_avb4.h"
#include "cluster/include/compute_add_avb.h"
#include "cluster/include/compute_remove_avb.h"
#include "cluster/include/trial_avb2.h"
#include "cluster/include/trial_avb4.h"
#include "gibbs/include/compute_gibbs_particle_transfer.h"
#include "chain/include/perturb_pivot.h"
#include "chain/include/select_segment.h"
#include "chain/include/perturb_crankshaft.h"
#include "chain/include/select_reptate.h"
#include "chain/include/perturb_reptate.h"
#include "chain/include/perturb_connector.h"
#include "chain/include/perturb_distance_angle_connector.h"
#include "chain/include/select_site_of_type.h"
#include "chain/include/perturb_site_type.h"
#include "chain/include/trial_grow.h"
#include "chain/include/select_branch.h"
#include "chain/include/perturb_branch.h"
#include "chain/include/perturb_to_anchor.h"
#include "chain/include/select_two_sites.h"
#include "chain/include/perturb_position_swap.h"

namespace feasst {

class MapTrialGrow {
 public:
  MapTrialGrow() {
    auto obj = MakeTrialGrow({{{"particle_type", "0"}, {"bond", "1"}, {"mobile_site", "1"}, {"anchor_site", "0"}}});
    obj->deserialize_map()["TrialGrow"] = obj;
  }
};

static MapTrialGrow mapper_ = MapTrialGrow();

void TrialGrow::build_(std::vector<argtype> * args) {
  const int num_args = static_cast<int>(args->size());
  ASSERT(num_args > 0, "TrialGrow requires args.size: " << num_args << " > 0");
  double weight = dble("weight", &(*args)[0], -1.);
  double weight_per_number_fraction = dble("weight_per_number_fraction", &(*args)[0], -1.);
  double number_fraction_exclude_type = integer("number_fraction_exclude_type", &(*args)[0], -1.);
  ASSERT(weight < 0 || weight_per_number_fraction < 0,
    "TrialGrow cannot have both weight and weight_per_number_fraction arguments.");
  const std::string default_print_num_accepted  = str("print_num_accepted", &(*args)[0], "false");
  const std::string default_num_steps = str("default_num_steps", &(*args)[0], "1");
  const std::string default_reference_index = str("default_reference_index", &(*args)[0], "-1");
  const std::string default_new_only = str("default_new_only", &(*args)[0], "false");
  // First, determine all trial types from args[0]
  std::vector<std::string> trial_types;
  std::vector<bool> trial_half_weight;
  if (used("add", (*args)[0])) {
    str("add", &(*args)[0]);
    trial_types.push_back("add");
    trial_half_weight.push_back(false);
  }
  if (used("remove", (*args)[0])) {
    str("remove", &(*args)[0]);
    trial_types.push_back("remove");
    trial_half_weight.push_back(false);
  }
  if (used("transfer", (*args)[0])) {
    str("transfer", &(*args)[0]);
    trial_types.push_back("add");
    trial_types.push_back("remove");
    trial_half_weight.push_back(true);
    trial_half_weight.push_back(true);
  }
  if (used("gibbs_transfer", (*args)[0])) {
    str("gibbs_transfer", &(*args)[0]);
    trial_types.push_back("gibbs_transfer");
    trial_half_weight.push_back(true);
    trial_half_weight.push_back(true);
  }
  if (used("add_avb", (*args)[0])) {
    str("add_avb", &(*args)[0]);
    trial_types.push_back("add_avb");
    trial_half_weight.push_back(false);
  }
  if (used("remove_avb", (*args)[0])) {
    str("remove_avb", &(*args)[0]);
    trial_types.push_back("remove_avb");
    trial_half_weight.push_back(false);
  }
  if (used("transfer_avb", (*args)[0])) {
    str("transfer_avb", &(*args)[0]);
    trial_types.push_back("add_avb");
    trial_types.push_back("remove_avb");
    trial_half_weight.push_back(true);
    trial_half_weight.push_back(true);
  }
  if (used("regrow_avb2", (*args)[0])) {
    str("regrow_avb2", &(*args)[0]);
    trial_types.push_back("regrow_avb2_in");
    trial_types.push_back("regrow_avb2_out");
    trial_half_weight.push_back(true);
    trial_half_weight.push_back(true);
  }
  if (used("regrow", (*args)[0])) {
    str("regrow", &(*args)[0]);
    trial_types.push_back("regrow");
    trial_half_weight.push_back(false);
  }
  if (used("regrow_avb4", (*args)[0])) {
    str("regrow_avb4", &(*args)[0]);
    trial_types.push_back("regrow_avb4");
    trial_half_weight.push_back(false);
  }
  if (used("translate", (*args)[0])) {
    str("translate", &(*args)[0]);
    trial_types.push_back("translate");
    trial_half_weight.push_back(false);
  }
  if (trial_types.size() == 0) {
    trial_types.push_back("partial_regrow");
    trial_half_weight.push_back(false);
  }

  // Second, determine the particle_type and first site from args[0]
  const std::string site = str("site", &(*args)[0], "-1");
  const std::string particle_type = str("particle_type", &(*args)[0]);

  // Third, determine the configuration_index from args[0]
  const std::string configuration_index = str("configuration_index", &(*args)[0], "0");
  const int configuration_index2 = integer("configuration_index2", &(*args)[0], -1);

  // Finally, add each trial to the factory
  for (int itr = 0; itr < static_cast<int>(trial_types.size()); ++itr) {
    const std::string& trial_type = trial_types[itr];
    DEBUG("trial_type: " << trial_type);
    std::shared_ptr<Trial> trial = MakeTrial({{"print_num_accepted",
                                        default_print_num_accepted}});
    trial->set_description("TrialGrow" + trial_type);
    if (weight > 0) {
      trial->set_weight(weight);
    }
    if (weight_per_number_fraction > 0) {
      trial->set_weight_per_number_fraction(weight_per_number_fraction);
    }
    if (number_fraction_exclude_type > 0) {
      trial->set_number_fraction_exclude_type(number_fraction_exclude_type);
    }
    if (trial_half_weight[itr]) {
      if (weight_per_number_fraction > 0) {
        trial->set_weight_per_number_fraction(trial->weight_per_number_fraction()/2.);
      } else {
        trial->set_weight(trial->weight()/2.);
      }
    }
    std::shared_ptr<TrialCompute> compute;
    for (int iarg = 0; iarg < num_args; ++iarg) {
      DEBUG("iarg: " << iarg);
      argtype iargs = (*args)[iarg];
      std::shared_ptr<TrialSelect> select;
      std::shared_ptr<Perturb> perturb;
      if (iarg == 0 && trial_type != "partial_regrow") {
        argtype sel_args = {{"particle_type", particle_type},
                            {"configuration_index", configuration_index},
                            {"site", site},
                            {"max_particles", str("max_particles", &iargs, "-1")},
                            {"min_particles", str("min_particles", &iargs, "-1")}};
        select = MakeTrialSelectParticle(sel_args);
        if (trial_type == "translate") {
          perturb = std::make_shared<PerturbTranslate>(&iargs);
          compute = std::make_shared<TrialComputeTranslate>(&iargs);
        } else if (trial_type == "add") {
          perturb = MakePerturbAdd();
          compute = MakeTrialComputeAdd();
        } else if (trial_type == "remove") {
          perturb = MakePerturbRemove();
          compute = MakeTrialComputeRemove();
        } else if (trial_type == "gibbs_transfer") {
          perturb = MakePerturbAdd();
          compute = MakeComputeGibbsParticleTransfer();
        } else if (trial_type == "regrow") {
          perturb = MakePerturbAnywhere();
          compute = MakeTrialComputeMove();
        } else if (trial_type == "add_avb") {
          iargs.insert({"particle_type", particle_type});
          iargs.insert({"configuration_index", configuration_index});
          iargs.insert({"site", site});
          iargs.insert({"grand_canonical", "true"});
          select = std::make_shared<SelectParticleAVB>(&iargs);
          auto add_avb = std::make_shared<PerturbAddAVB>(&iargs);
          ASSERT(add_avb->delay_add(), "ComputeAddAVB assumes delay_add");
          perturb = add_avb;
          compute = std::make_shared<ComputeAddAVB>();
        } else if (trial_type == "remove_avb") {
          iargs.insert({"particle_type", particle_type});
          iargs.insert({"configuration_index", configuration_index});
          iargs.insert({"site", site});
          iargs.insert({"grand_canonical", "true"});
          select = std::make_shared<SelectParticleAVB>(&iargs);
          perturb = std::make_shared<PerturbRemove>(
            MakePerturbMoveAVB({{{"inside", "true"}}}));
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
          iargs.insert({"configuration_index", configuration_index});
          iargs.insert({"site", site});
          select = std::make_shared<SelectParticleAVB>(&iargs);
          perturb = std::make_shared<PerturbMoveAVB>(&perturb_args);
          compute = std::make_shared<ComputeAVB2>(&iargs);
        } else if (trial_type == "regrow_avb4") {
          //argtype args_sel, args_mv;
          gen_avb4_args_(&iargs);
          iargs.insert({"particle_type", particle_type});
          iargs.insert({"configuration_index", configuration_index});
          iargs.insert({"site", site});
          select = std::make_shared<SelectParticleAVB>(&iargs);
          perturb = std::make_shared<PerturbMoveAVB>(&iargs);
          compute = std::make_shared<ComputeAVB4>();
        } else {
          FATAL("unreocgnized trial_type: " << trial_type);
        }
      } else {
        int used = 0;
        if (boolean("bond", &iargs, false)) {
          ++used;
          select = MakeTrialSelectBond({
            {"particle_type", particle_type},
            {"configuration_index", configuration_index},
            {"mobile_site", str("mobile_site", &iargs)},
            {"anchor_site", str("anchor_site", &iargs)}});
          perturb = std::make_shared<PerturbDistance>(&iargs);
        }
        if (boolean("angle", &iargs, false)) {
          ASSERT(used == 0, "cannot have more than one");
          ++used;
          select = MakeTrialSelectAngle({
            {"particle_type", particle_type},
            {"configuration_index", configuration_index},
            {"mobile_site", str("mobile_site", &iargs)},
            {"anchor_site", str("anchor_site", &iargs)},
            {"anchor_site2", str("anchor_site2", &iargs)}});
          perturb = std::make_shared<PerturbDistanceAngle>(&iargs);
        }
        if (boolean("dihedral", &iargs, false)) {
          ASSERT(used == 0, "cannot have more than one");
          ++used;
          select = MakeTrialSelectDihedral({
            {"particle_type", particle_type},
            {"configuration_index", configuration_index},
            {"mobile_site", str("mobile_site", &iargs)},
            {"anchor_site", str("anchor_site", &iargs)},
            {"anchor_site2", str("anchor_site2", &iargs)},
            {"anchor_site3", str("anchor_site3", &iargs)}});
          perturb = std::make_shared<PerturbDihedral>(&iargs);
        }
        if (boolean("branch", &iargs, false)) {
          ASSERT(used == 0, "cannot have more than one");
          ++used;
          select = MakeSelectBranch({
            {"particle_type", particle_type},
            {"configuration_index", configuration_index},
            {"mobile_site", str("mobile_site", &iargs)},
            {"mobile_site2", str("mobile_site2", &iargs)},
            {"anchor_site", str("anchor_site", &iargs)},
            {"anchor_site2", str("anchor_site2", &iargs)}});
          perturb = std::make_shared<PerturbBranch>(&iargs);
        }
        if (boolean("reptate", &iargs, false)) {
          FATAL("TrialGrow::reptate has an issue and should not be used.");
          ASSERT(used == 0, "cannot have more than one");
          ++used;
          select = MakeTrialSelectBond({
            {"particle_type", particle_type},
            {"configuration_index", configuration_index},
            {"mobile_site", str("mobile_site", &iargs)},
            {"anchor_site", str("anchor_site", &iargs)},
            {"ignore_bond", "true"}});
          perturb = std::make_shared<PerturbToAnchor>(&iargs);
        }
        if (boolean("position_swap", &iargs, false)) {
          ASSERT(used == 0, "cannot have more than one");
          ++used;
          select = MakeSelectTwoSites({
            {"particle_type", particle_type},
            {"particle_type2", str("particle_type2", &iargs)},
            {"configuration_index", configuration_index},
            {"mobile_site", str("mobile_site", &iargs)},
            {"mobile_site2", str("mobile_site2", &iargs)}});
          perturb = std::make_shared<PerturbPositionSwap>(&iargs);
        }
        if (boolean("rigid_body_connector", &iargs, false)) {
          ASSERT(used == 0, "cannot have more than one");
          ++used;
          select = MakeTrialSelectBond({
            {"particle_type", particle_type},
            {"configuration_index", configuration_index},
            {"mobile_site", str("mobile_site", &iargs)},
            {"anchor_site", str("anchor_site", &iargs)}});
          perturb = std::make_shared<PerturbConnector>(&iargs);
        }
        if (boolean("rigid_body_angle", &iargs, false)) {
          ASSERT(used == 0, "cannot have more than one");
          ++used;
          select = MakeTrialSelectAngle({
            {"particle_type", particle_type},
            {"configuration_index", configuration_index},
            {"mobile_site", str("mobile_site", &iargs)},
            {"anchor_site", str("anchor_site", &iargs)},
            {"anchor_site2", str("anchor_site2", &iargs)}});
          perturb = std::make_shared<PerturbDistanceAngleConnector>(&iargs);
        }
        ASSERT(used == 1, "args: " << str(iargs) <<
          ". Requires one of bond, angle, dihedral, branch, reptate, etc");
        if (!compute) {
          DEBUG("num_args " << num_args);
          // HWH: disable tunable trials until developed properly
//          if (num_args == 1) {
//            compute = std::make_shared<TrialComputeTranslate>(&iargs);
//          } else {
//            compute = std::make_shared<TrialComputeMove>(&iargs);
//          }
          compute = std::make_shared<TrialComputeMove>(&iargs);
        } else {
          DEBUG(compute->class_name());
        }
      }
      const std::string num_steps = str("num_steps", &iargs, default_num_steps);
      //ASSERT(trial_type != "translate" || num_steps == "1",
      //  "For " << trial_type << ", num_steps must be 1");
      argtype stage_args = {{"num_steps", num_steps},
        {"reference_index", str("reference_index", &iargs, default_reference_index)},
        {"new_only", str("new_only", &iargs, default_new_only)},
      };
      feasst_check_all_used(iargs);
      trial->add_stage(select, perturb, &stage_args);
      feasst_check_all_used(stage_args);
    }
    if (trial_types[0] == "gibbs_transfer") {
      ASSERT(static_cast<int>(trial_types.size()) == 1,
        "gibbs_transfer should only have one trial.");
      ASSERT(configuration_index2 > 0,
        "configuration_index2 is a required argument for gibbs_transfer.");
      // duplicate all stages again except, for the first new stage, replace
      // PerturbAdd with PerturbRemove and select the other configuration.
      const int num_stages = trial->num_stages();
      for (int istage = 0; istage < num_stages; ++istage) {
        DEBUG("istage " << istage);
        const TrialStage& stage = trial->stage(istage);
        DEBUG("pert name " << stage.perturb().class_name());
        auto cp_stg = std::make_shared<TrialStage>(deep_copy(stage));
        cp_stg->get_trial_select()->set_configuration_index(configuration_index2);
        trial->add_stage(cp_stg);
        if (istage == 0) {
          TrialStage * stg0 = trial->get_stage_(num_stages);
          TrialSelect * sel0 = stg0->get_trial_select();
          sel0->set_ghost(false);
          stg0->set(MakePerturbRemove());
        }
      }
//      for (std::shared_ptr<TrialStage> stage : trial->stages()) {
//        DEBUG(stage->perturb().class_name());
//      }
    }
    trial->set(compute);
    DEBUG("compute " << compute->class_name());
    add(trial);
    if (trial_types[0] == "gibbs_transfer") {
      // duplicate gibbs_transfer with reverse configuration_index
      auto cp_trl = std::make_shared<Trial>(deep_copy(*trial));
      const int conf_ind1 = str_to_int(configuration_index);
      const int conf_ind2 = configuration_index2;
      for (int istage = 0; istage < cp_trl->num_stages(); ++istage) {
//        TrialStage * stage = cp_trl->get_stage_(istage);
//        TrialSelect * sel = stage->get_trial_select();
        TrialSelect * sel = cp_trl->get_stage_(istage)->get_trial_select();
        DEBUG("istage " << istage << " conf " << sel->configuration_index());
        if (sel->configuration_index() == conf_ind1) {
          sel->set_configuration_index(conf_ind2);
        } else if (sel->configuration_index() == conf_ind2) {
          sel->set_configuration_index(conf_ind1);
        } else {
          FATAL("unrecognized conf index: " << sel->configuration_index());
        }
      }
      add(cp_trl);
    }
  }
}

TrialGrow::TrialGrow(std::vector<argtype> args) : TrialFactoryNamed() {
  class_name_ = "TrialGrow";
  build_(&args);
}

class MapTrialGrowFile {
 public:
  MapTrialGrowFile() {
    auto obj = std::make_shared<TrialGrowFile>();
    obj->deserialize_map()["TrialGrowFile"] = obj;
  }
};

static MapTrialGrowFile mapper_trial_grow_file_ = MapTrialGrowFile();

void TrialGrowFile::add_(const argtype add_args, std::vector<argtype> * args) {
  argtype * arg0 = &(*args)[0];
  arg0->insert(add_args.begin(), add_args.end());
  build_(args);
  args->clear();
}

// Convert file contents to a list of argtypes to use in the above constructor
TrialGrowFile::TrialGrowFile(argtype * args) : TrialGrow() {
  class_name_ = "TrialGrowFile";
  std::string grow_file = str("grow_file", args, "");
  if (grow_file.empty()) {
    if (used("file_name", *args)) {
      WARN("TrialGrowFile::file_name renamed to grow_file.");
      grow_file = str("file_name", args);
    } else {
      FATAL("TrialGrowFile::grow_file is a required argument.");
    }
  }
  std::vector<argtype> reformated;
  std::ifstream file(grow_file);
  ASSERT(file.good(), "cannot find " << grow_file);
  // until end of file, find TrialGrowFile
  // check for optional empty line
  // read each line as all arguments in one stage
  // when empty line or end of file is reached, add Trial and repeat
  const bool is_found = find("TrialGrowFile", file);
  ASSERT(is_found, "TrialGrowFile not found in " << grow_file);
  std::string line;
  while (std::getline(file, line)) {
    if (line.empty()) {
      if (reformated.size() > 0) {
        add_(*args, &reformated);
      }
    } else {
      if (line[0] != '#') {
        reformated.push_back(line_to_argtype(line));
      }
    }
  }
  if (reformated.size() > 0) {
    add_(*args, &reformated);
  }
}
TrialGrowFile::TrialGrowFile(argtype args) : TrialGrowFile(&args) {
  feasst_check_all_used(args);
}

}  // namespace feasst
