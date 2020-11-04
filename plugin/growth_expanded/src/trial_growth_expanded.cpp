#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/random.h"
#include "configuration/include/domain.h"
#include "monte_carlo/include/trials.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "growth_expanded/include/trial_growth_expanded.h"

namespace feasst {

class MapTrialComputeGrowAdd {
 public:
  MapTrialComputeGrowAdd() {
    auto obj = std::make_shared<TrialComputeGrowAdd>();
    obj->deserialize_map()["TrialComputeGrowAdd"] = obj;
  }
};

static MapTrialComputeGrowAdd mapper_trial_compute_grow_add_ = MapTrialComputeGrowAdd();

std::shared_ptr<TrialCompute> TrialComputeGrowAdd::create(std::istream& istr) const {
  return std::make_shared<TrialComputeGrowAdd>(istr);
}

TrialComputeGrowAdd::TrialComputeGrowAdd(std::istream& istr)
  : TrialCompute(istr) {
  ASSERT(class_name_ == "TrialComputeGrowAdd", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(154 == version, "mismatch version: " << version);
}


void TrialComputeGrowAdd::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_compute_(ostr);
  feasst_serialize_version(154, ostr);
}

class MapTrialComputeGrowRemove {
 public:
  MapTrialComputeGrowRemove() {
    auto obj = std::make_shared<TrialComputeGrowRemove>();
    obj->deserialize_map()["TrialComputeGrowRemove"] = obj;
  }
};

static MapTrialComputeGrowRemove mapper_trial_compute_grow_remove_ = MapTrialComputeGrowRemove();

std::shared_ptr<TrialCompute> TrialComputeGrowRemove::create(std::istream& istr) const {
  return std::make_shared<TrialComputeGrowRemove>(istr);
}

TrialComputeGrowRemove::TrialComputeGrowRemove(std::istream& istr)
  : TrialCompute(istr) {
  ASSERT(class_name_ == "TrialComputeGrowRemove", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(225 == version, "mismatch version: " << version);
}


void TrialComputeGrowRemove::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_compute_(ostr);
  feasst_serialize_version(225, ostr);
}

class MapTrialComputeGrow {
 public:
  MapTrialComputeGrow() {
    auto obj = std::make_shared<TrialComputeGrow>();
    obj->deserialize_map()["TrialComputeGrow"] = obj;
  }
};

static MapTrialComputeGrow mapper_trial_compute_grow_ = MapTrialComputeGrow();

std::shared_ptr<TrialCompute> TrialComputeGrow::create(std::istream& istr) const {
  return std::make_shared<TrialComputeGrow>(istr);
}

TrialComputeGrow::TrialComputeGrow(std::istream& istr)
  : TrialCompute(istr) {
  ASSERT(class_name_ == "TrialComputeGrow", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(155 == version, "mismatch version: " << version);
}


void TrialComputeGrow::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_compute_(ostr);
  feasst_serialize_version(155, ostr);
}

class MapTrialGrowthExpanded {
 public:
  MapTrialGrowthExpanded() {
    auto obj = std::make_shared<TrialGrowthExpanded>(MakeTrialAdd(
      {{"particle_type", "0"}}), MakeTrialAdd({{"particle_type", "0"}}));
    obj->deserialize_map()["TrialGrowthExpanded"] = obj;
  }
};

static MapTrialGrowthExpanded mapper_trial_growth_expanded_ = MapTrialGrowthExpanded();

std::shared_ptr<Trial> TrialGrowthExpanded::create(std::istream& istr) const {
  return std::make_shared<TrialGrowthExpanded>(istr);
}

TrialGrowthExpanded::TrialGrowthExpanded(std::istream& istr)
  : Trial(istr) {
  ASSERT(class_name_ == "TrialGrowthExpanded", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(793 == version, "mismatch version: " << version);
}


void TrialGrowthExpanded::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(793, ostr);
}

TrialGrowthExpanded::TrialGrowthExpanded(
    std::shared_ptr<Trial> grow,
    std::shared_ptr<Trial> shrink,
    const argtype& args) : Trial(args) {
  class_name_ = "TrialGrowthExpanded";
  grow_ = grow;
  shrink_ = shrink;
  compute_add_ = std::make_shared<TrialComputeGrowAdd>();
  compute_remove_ = std::make_shared<TrialComputeGrowRemove>();
  compute_grow_ = std::make_shared<TrialComputeGrow>();
  add_(grow_->stages()[0]);
  DEBUG(shrink_);
  DEBUG("num " << shrink_->stages().size());
  const TrialSelect& sel = shrink_->stages()[0]->trial_select();
  growing_particle_ = MakeTrialSelectParticle({
    {"particle_type", str(sel.particle_type())},
    {"site", "0"}, // HWH hardcoded for site0
  });

  // the first selection stage of shrink should use perturbed
  auto stage = shrink_->stages()[0];
  stage->set(std::make_shared<SelectPerturbed>());
  shrink_->set(0, stage);
}

void TrialGrowthExpanded::before_select(Acceptance * acceptance, Criteria * criteria) {
  criteria->set_expanded_state(current_growth_stage_(growing_),
                            num_growth_stages());
  if (growing_) {
    DEBUG("attempt grow, stage: " << growth_stage_);
    set(0, grow_->stages()[growth_stage_]);
    if (growth_stage_ == 0) {
      set(compute_add_);
    } else {
      compute_grow_->set_grow();
      set(compute_grow_);
    }
  } else {
    DEBUG("attempt shrink, stage: " << growth_stage_);
    if (growth_stage_ == 0) {
      set(0, shrink_->stages()[num_growth_stages() - 1]);
      compute_grow_->set_shrink();
      set(compute_grow_);
    } else if (growth_stage_ == 1) {
      set(0, shrink_->stages()[0]);
      set(compute_remove_);
    } else {
      set(0, shrink_->stages()[growth_stage_ - 1]);
      compute_grow_->set_shrink();
      set(compute_grow_);
    }
  }

  // acceptanced::perturbed is used for selection of bonds, etc.
  DEBUG("growingp " << growing_particle_->mobile().str());
  if (growth_stage_ != 0) {
    acceptance->add_to_perturbed(growing_particle_->mobile());
  }
}

void TrialComputeGrowAdd::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("TrialComputeGrowAdd");
  compute_rosenbluth(0, criteria, system, acceptance, stages, random);
  const TrialSelect& select = (*stages)[0]->trial_select();
  system->get_configuration()->revive(select.mobile());
  acceptance->set_energy_new(criteria->current_energy() + acceptance->energy_new());
  { // Metropolis
    const Configuration& config = system->configuration();
    const double volume = config.domain().volume();
    const int particle_index = select.mobile().particle_index(0);
    const int particle_type = config.select_particle(particle_index).type();
    const ThermoParams& params = system->thermo_params();
    DEBUG("volume " << volume << " selprob " << select.probability() <<
      " betamu " << params.beta_mu(particle_type));
    const double power = 1.;///static_cast<double>(criteria->num_expanded_states());
    acceptance->add_to_ln_metropolis_prob(
      power*(std::log(volume)
             + params.beta_mu(particle_type)*criteria->num_expanded_states())
    );
  }
}

void TrialComputeGrowRemove::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("TrialComputeRemove");
  compute_rosenbluth(1, criteria, system, acceptance, stages, random);
  acceptance->set_energy_new(criteria->current_energy() - acceptance->energy_old());
  acceptance->add_to_macrostate_shift(-1);
  { // Metropolis
    const Configuration& config = system->configuration();
    const double volume = config.domain().volume();
    const TrialSelect& select = (*stages)[0]->trial_select();
    const int particle_index = select.mobile().particle_index(0);
    const int particle_type = config.select_particle(particle_index).type();
    const ThermoParams& params = system->thermo_params();
    DEBUG("volume " << volume << " selprob " << select.probability() << " betamu " << params.beta_mu(particle_type));
    const double power = 1;//./static_cast<double>(criteria->num_expanded_states());
    acceptance->add_to_ln_metropolis_prob(
      power*(-std::log(volume) - params.beta_mu(particle_type)));
    DEBUG("lnmet " << acceptance->ln_metropolis_prob());
  }
}

bool TrialGrowthExpanded::attempt(Criteria * criteria,
    System * system,
    Random * random) {
  DEBUG("num " << system->configuration().num_particles());
  // whether accepted or rejected, select new growing particle when stage0
  if (growth_stage_ == 0) {
    growing_particle_->select(growing_particle_->anchor(), system, random);
  }
  growing_ = random->coin_flip();
  const bool accepted = Trial::attempt(criteria, system, random);
  DEBUG("accepted? " << accepted);
  if (accepted) {
    growth_stage_ = current_growth_stage_(growing_);
    update_growing_particle_();
  }
  if ( (accepted and !growing_) or (!accepted and growing_) ) {
    get_stage_(0)->set_mobile_physical(false, system);
  }
  criteria->set_expanded_state(growth_stage_, num_growth_stages());
  DEBUG("growingpend " << growing_particle_->mobile().str());
  return accepted;
}

}  // namespace feasst
