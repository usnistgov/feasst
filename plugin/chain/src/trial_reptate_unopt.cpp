#include "utils/include/debug.h"
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/position.h"
#include "configuration/include/particle.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "monte_carlo/include/perturb_dihedral.h"
#include "monte_carlo/include/perturb_distance_angle.h"
#include "monte_carlo/include/perturb_distance.h"
#include "monte_carlo/include/trial_select_dihedral.h"
#include "monte_carlo/include/trial_select_angle.h"
#include "monte_carlo/include/trial_select_bond.h"
#include "monte_carlo/include/trial_compute_move.h"
#include "chain/include/perturb_to_anchor.h"
#include "chain/include/trial_reptate_unopt.h"

namespace feasst {

FEASST_MAPPER(TrialReptateUnoptHalf,);

TrialReptateUnoptHalf::TrialReptateUnoptHalf(argtype * args) : Trial(args) {
  class_name_ = "TrialReptateUnoptHalf";
  // store and process args at precompute
  args_ = *args;
  args->clear();
}
TrialReptateUnoptHalf::TrialReptateUnoptHalf(argtype args) : TrialReptateUnoptHalf(&args) {
  feasst_check_all_used(args);
}
TrialReptateUnoptHalf::~TrialReptateUnoptHalf() {}

void TrialReptateUnoptHalf::precompute(Criteria * criteria, System * system) {
  FATAL("This trial is disabled because PerturbToAnchor does not properly "
    << "account for excluded bond energies. In addition, Perturb needs to "
    << "act differently for old vs new configs (e.g., PerturbToAnchor in new "
    << "config and Perturb[Bond/Angle/Dihedral] in the old one, especially "
    << "if num_steps > 1 is enabled. Bond excluded energies in the Rosenbluth "
    << "may need to be refactored to make this easier. Perturb may also need "
    << "to be refactored to allow different ones.");
  DEBUG("Parse args");
  const int ptype = integer("particle_type", &args_);
  const int config_index = integer("configuration_index", &args_, 0);
  const std::string reference_index = str("reference_index", &args_, "-1");
  ASSERT(feasst::str_to_int(reference_index) > -1, "TrialReptateUnopt requires a reference_index");
  const bool reverse = boolean("reverse", &args_);
  const int num_steps = integer("num_steps", &args_, 1);
  ASSERT(num_steps == 1, "TrialReptateUnoptHalf requires num_steps == 1");
  feasst_check_all_used(args_);

  const Configuration& config = system->configuration(config_index);
  const Particle& part = config.particle_type(ptype);
  const int num_sites = part.num_sites();
  DEBUG("num_sites " << part.num_sites());
  argtype sel_begin_args = {{"particle_type", str(ptype)}, {"ignore_bond", "true"},
    {"configuration_index", str(config_index)}};
  argtype sel_end_args = {{"particle_type", str(ptype)}};
  if (reverse) {
    sel_begin_args.insert({"mobile_site", str(num_sites-1)});
    sel_begin_args.insert({"anchor_site", str(num_sites-2)});
    sel_end_args.insert({"mobile_site", "0"});
    sel_end_args.insert({"anchor_site", "1"});
  } else {
    sel_begin_args.insert({"mobile_site", "0"});
    sel_begin_args.insert({"anchor_site", "1"});
    sel_end_args.insert({"mobile_site", str(num_sites-1)});
    sel_end_args.insert({"anchor_site", str(num_sites-2)});
  }
  auto sel_begin = std::make_shared<TrialSelectBond>(sel_begin_args);
  auto pert_begin = std::make_shared<PerturbToAnchor>();
  if (part.num_angles() > 0) {
    if (reverse) {
      sel_end_args.insert({"anchor_site2", "2"});
    } else {
      sel_end_args.insert({"anchor_site2", str(num_sites-3)});
    }
  }
  if (part.num_dihedrals() > 0) {
    if (reverse) {
      sel_end_args.insert({"anchor_site3", "3"});
    } else {
      sel_end_args.insert({"anchor_site3", str(num_sites-4)});
    }
  }
  std::shared_ptr<Perturb> pert_end;
  std::shared_ptr<TrialSelect> sel_end;
  if (part.num_dihedrals() > 0) {
    sel_end = std::make_shared<TrialSelectDihedral>(sel_end_args);
    pert_end = std::make_shared<PerturbDihedral>();
  } else if (part.num_angles() > 0) {
    sel_end = std::make_shared<TrialSelectAngle>(sel_end_args);
    pert_end = std::make_shared<PerturbDistanceAngle>();
  } else {
    sel_end = std::make_shared<TrialSelectBond>(sel_end_args);
    pert_end = std::make_shared<PerturbDistance>();
  }
  argtype stage_args = {{"reference_index", reference_index}};
  add_stage(sel_begin, pert_begin, &stage_args);
  if (reverse) {
    for (int site = num_sites - 2; site > 0; --site) {
      stage_args = {{"reference_index", reference_index}};
      add_stage(
        std::make_shared<TrialSelectBond>(argtype({
          {"particle_type", str(ptype)},
          {"mobile_site", str(site)},
          {"anchor_site", str(site - 1)},
          {"ignore_bond", "true"},
        })),
        std::make_shared<PerturbToAnchor>(),
        &stage_args
      );
    }
  } else {
    for (int site = 1; site < num_sites - 1; ++site) {
      add_stage(
        std::make_shared<TrialSelectBond>(argtype({
          {"particle_type", str(ptype)},
          {"mobile_site", str(site)},
          {"anchor_site", str(site + 1)},
          {"ignore_bond", "true"},
        })),
        std::make_shared<PerturbToAnchor>(),
        &stage_args
      );
    }
  }
  stage_args = {{"reference_index", reference_index}};
  add_stage(sel_end, pert_end, &stage_args);
  set(std::make_shared<TrialComputeMove>());
  Trial::precompute(criteria, system);
}

TrialReptateUnoptHalf::TrialReptateUnoptHalf(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 1034, "mismatch version: " << version);
}

void TrialReptateUnoptHalf::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(1034, ostr);
}

FEASST_MAPPER(TrialReptateUnopt,);

TrialReptateUnopt::TrialReptateUnopt(argtype * args) : TrialFactoryNamed() {
  class_name_ = "TrialReptateUnopt";
  argtype reverse_args(*args);
  reverse_args.insert({"reverse", "true"});
  args->insert({"reverse", "false"});
  auto reverse_trial = std::make_shared<TrialReptateUnoptHalf>(reverse_args);
  reverse_trial->set_weight(reverse_trial->weight()/2.);
  add(reverse_trial);
  auto trial = std::make_shared<TrialReptateUnoptHalf>(args);
  trial->set_weight(trial->weight()/2.);
  add(trial);
}
TrialReptateUnopt::TrialReptateUnopt(argtype args) : TrialReptateUnopt(&args) {
  feasst_check_all_used(args);
}

}  // namespace feasst
