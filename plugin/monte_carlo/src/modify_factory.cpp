#include "utils/include/serialize.h"
#include "monte_carlo/include/modify_factory.h"

namespace feasst {

class MapModifyFactory {
 public:
  MapModifyFactory() {
    ModifyFactory().deserialize_map()["ModifyFactory"] =
      std::make_shared<ModifyFactory>();
  }
};

static MapModifyFactory mapper_ = MapModifyFactory();

void ModifyFactory::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  for (const std::shared_ptr<Modify>& modify : modifiers_) {
    modify->initialize(criteria, system, trial_factory);
  }
}

void ModifyFactory::trial_(Criteria* criteria,
    System* system,
    Random * random,
    TrialFactory* trial_factory,
    const int index) {
  // timer_.start(index + 1);
  DEBUG("index " << index << " sz " << modifiers_.size());
  ASSERT(index < static_cast<int>(modifiers_.size()),
    "index: " << index << " too large when there are " << modifiers_.size());
  DEBUG(modifiers_[index]->class_name());
  modifiers_[index]->trial(criteria, system, random, trial_factory);
  // timer_.end();
}

void ModifyFactory::trial(Criteria * criteria,
    System * system,
    Random * random,
    TrialFactory * trial_factory) {
  DEBUG(" class? " << class_name());
  int stt = -1;
  if (is_multistate()) {
    stt = criteria->state();
  }
  if ( (stop_after_phase() != -1 && criteria->phase() > stop_after_phase()) ||
       (stop_after_iteration() != -1 && criteria->num_iterations(stt) > stop_after_iteration()) ||
       (criteria->phase() <= start_after_phase()) ||
       (criteria->num_iterations(stt) <= start_after_iteration()) ) {
    return;
  }
  if (is_multistate()) {
    DEBUG("multistate");
    DEBUG("state? " << criteria->state());
    if (is_multistate_aggregate()) {
      DEBUG("aggregating");
      DEBUG("sz " << modifiers_.size());
      ASSERT(criteria->state() < static_cast<int>(modifiers_.size()),
        "state: " << criteria->state() << " >= multistate modifiers: " <<
        modifiers_.size() << ". Was a flat histogram simulation reinitialized"
        << " after a multistate Modifier?");
      DEBUG(modifiers_[criteria->state()]->class_name());
      modifiers_[criteria->state()]->check_update_(criteria, system, random, trial_factory);
      DEBUG("is time? " << trials_per_write() << " " << trials_since_write_);
      if (is_time(trials_per_write(), &trials_since_write_)) {
        std::stringstream ss;
        for (int state = 0; state < num(); ++state) {
          if (state == 0) {
            ss << "state,";
            ss << modifiers_[state]->header(*criteria, *system, *trial_factory);
          }
          DEBUG("state " << state);
          DEBUG("crit " << criteria->state());
          DEBUG("crit " << criteria->num_states());
          ss << state << ",";
          ss << modifiers_[state]->write(criteria, system, trial_factory);
//          ss << std::endl;
        }
        printer(ss.str(), file_name(*criteria));
      }
    } else {
      trial_(criteria, system, random, trial_factory, criteria->state());
    }
  } else {
    DEBUG("not multistate");
    for (int index = 0; index < num(); ++index) {
      DEBUG("index " << index);
      trial_(criteria, system, random, trial_factory, index);
    }
  }
}

void ModifyFactory::write_to_file(Criteria * criteria,
  System * system,
  TrialFactory * trial_factory) {
  if (is_multistate()) {
    if (is_multistate_aggregate()) {
      std::stringstream ss;
      for (int state = 0; state < num(); ++state) {
        if (state == 0) {
          ss << "state,";
          ss << modifiers_[state]->header(*criteria, *system, *trial_factory);
        }
        ss << state << ",";
        ss << modifiers_[state]->write(criteria, system, trial_factory);
      }
      printer(ss.str(), file_name(*criteria));
    } else {
      modifiers_[criteria->state()]->write_to_file(criteria, system, trial_factory);
    }
  } else {
    for (int index = 0; index < num(); ++index) {
      modifiers_[index]->write_to_file(criteria, system, trial_factory);
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
}

void ModifyFactory::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(7177, ostr);
  feasst_serialize_fstdr(modifiers_, ostr);
}

}  // namespace feasst
