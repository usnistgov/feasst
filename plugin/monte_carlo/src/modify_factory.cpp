#include "utils/include/serialize.h"
#include "utils/include/timer_rdtsc.h"
#include "monte_carlo/include/trial_factory.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/modify_factory.h"

namespace feasst {

FEASST_MAPPER(ModifyFactory,);

ModifyFactory::ModifyFactory(argtype args) : Modify(&args) {}
ModifyFactory::~ModifyFactory() {}

void ModifyFactory::initialize(MonteCarlo * mc) {
  Modify::initialize(mc);
  for (const std::shared_ptr<Modify>& modify : modifiers_) {
    modify->initialize(mc);
  }
}

void ModifyFactory::add(std::shared_ptr<Modify> modify) {
  modifiers_.push_back(modify);
  if (timer_) timer_->add();
}

void ModifyFactory::remove(const int index) {
  modifiers_.erase(modifiers_.begin() + index);
  if (timer_) timer_->erase(index);
}

void ModifyFactory::trial_(MonteCarlo * mc, const int index) {
  DEBUG("index " << index << " sz " << modifiers_.size());
  ASSERT(index < static_cast<int>(modifiers_.size()),
    "index: " << index << " too large when there are " << modifiers_.size());
  DEBUG(modifiers_[index]->class_name());
  if (timer_) timer_->start(index);
  modifiers_[index]->trial(mc);
  if (timer_) timer_->start(-1);
}

void ModifyFactory::trial(MonteCarlo * mc) {
  DEBUG(" class? " << class_name());
  const Criteria& criteria = mc->criteria();
  int stt = -1;
  if (is_multistate()) {
    stt = criteria.state();
  }
  if ( (stop_after_phase() != -1 && criteria.phase() > stop_after_phase()) ||
       (stop_after_cycle() != -1 && criteria.num_cycles(stt) > stop_after_cycle()) ||
       (criteria.phase() <= start_after_phase()) ||
       (criteria.num_cycles(stt) <= start_after_cycle()) ) {
    return;
  }
  if (is_multistate()) {
    DEBUG("multistate");
    DEBUG("state? " << criteria.state());
    if (is_multistate_aggregate()) {
      DEBUG("aggregating");
      DEBUG("sz " << modifiers_.size());
      ASSERT(criteria.state() < static_cast<int>(modifiers_.size()),
        "state: " << criteria.state() << " >= multistate modifiers: " <<
        modifiers_.size() << ". Was a flat histogram simulation reinitialized"
        << " after a multistate Modifier?");
      DEBUG(modifiers_[criteria.state()]->class_name());
      modifiers_[criteria.state()]->check_update_(mc);
      DEBUG("is time? " << trials_per_write() << " " << trials_since_write_);
      if (is_time(trials_per_write(), &trials_since_write_)) {
        std::stringstream ss;
        for (int state = 0; state < num(); ++state) {
          if (state == 0) {
            ss << "state,";
            ss << modifiers_[state]->header(*mc);
          }
          DEBUG("state " << state);
          DEBUG("crit " << criteria.state());
          DEBUG("crit " << criteria.num_states());
          ss << state << ",";
          ss << modifiers_[state]->write(mc);
//          ss << std::endl;
        }
        printer(ss.str(), output_file(criteria));
      }
    } else {
      trial_(mc, criteria.state());
    }
  } else {
    DEBUG("not multistate");
    for (int index = 0; index < num(); ++index) {
      DEBUG("index " << index);
      trial_(mc, index);
    }
  }
}

void ModifyFactory::write_to_file(MonteCarlo * mc) {
  if (is_multistate()) {
    const Criteria& criteria = mc->criteria();
    if (is_multistate_aggregate()) {
      std::stringstream ss;
      for (int state = 0; state < num(); ++state) {
        if (state == 0) {
          ss << "state,";
          ss << modifiers_[state]->header(*mc);
        }
        ss << state << ",";
        ss << modifiers_[state]->write(mc);
      }
      printer(ss.str(), output_file(criteria));
    } else {
      modifiers_[criteria.state()]->write_to_file(mc);
    }
  } else {
    for (int index = 0; index < num(); ++index) {
      modifiers_[index]->write_to_file(mc);
    }
  }
}

void ModifyFactory::adjust_bounds(const bool adjusted_up,
    const std::vector<int>& states,
    ModifyFactory * modify_factory) {
  for (int ai = 0; ai < num(); ++ai) {
    if (modify(ai).is_multistate()) {
      DEBUG("ai " << ai);
      for (const int state : states) {
        if (adjusted_up) {
          *modify_factory->get_modify(ai)->get_modify(state) = modify(ai).modify(state);
        } else {
          *get_modify(ai)->get_modify(state) = modify_factory->modify(ai).modify(state);
        }
      }
    }
  }
}

void ModifyFactory::set_timer() {
  timer_ = std::make_unique<TimerRDTSC>(num());
}

void ModifyFactory::synchronize_(const Modify& modify) {
  Modify::synchronize_(modify);
  for (int imod = 0; imod < num(); ++imod) {
    modifiers_[imod]->synchronize_(modify.modify(imod));
  }
}

ModifyFactory::ModifyFactory(std::istream& istr) : Modify(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 7177, "unrecognized verison: " << version);
  // feasst_deserialize_fstdr(&modifiers_, istr);
  // HWH for unknown reasons, function template doesn't work
  int dim1;
  istr >> dim1;
  modifiers_.resize(dim1);
  for (int index = 0; index < dim1; ++index) {
    int existing;
    istr >> existing;
    if (existing != 0) {
      modifiers_[index] = modifiers_[index]->deserialize(istr);
    }
  }
  feasst_deserialize(timer_, istr);
}

void ModifyFactory::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(7177, ostr);
  feasst_serialize_fstdr(modifiers_, ostr);
  feasst_serialize(timer_, ostr);
}

}  // namespace feasst
