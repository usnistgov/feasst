
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "utils/include/timer.h"
#include "math/include/utils_math.h"
#include "system/include/visit_model_cell.h"
#include "monte_carlo/include/run.h"
#include "monte_carlo/include/monte_carlo.h"

namespace feasst {

Run::Run(argtype * args) {
  num_trials_ = integer("num_trials", args, -1);
  until_num_particles_ = integer("until_num_particles", args, -1);
  for_hours_ = dble("for_hours", args, -1);
  until_criteria_complete_ = boolean("until_criteria_complete", args, false);
  class_name_ = "Run";
}
Run::Run(argtype args) : Run(&args) {
  check_all_used(args);
}

class MapRun {
 public:
  MapRun() {
    auto obj = MakeRun();
    obj->deserialize_map()["Run"] = obj;
  }
};

static MapRun mapper_Run = MapRun();

Run::Run(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 3854, "mismatch version: " << version);
  feasst_deserialize(&num_trials_, istr);
  feasst_deserialize(&until_num_particles_, istr);
  feasst_deserialize(&for_hours_, istr);
  feasst_deserialize(&until_criteria_complete_, istr);
}

void Run::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(3854, ostr);
  feasst_serialize(num_trials_, ostr);
  feasst_serialize(until_num_particles_, ostr);
  feasst_serialize(for_hours_, ostr);
  feasst_serialize(until_criteria_complete_, ostr);
}

void Run::run(MonteCarlo * mc) {
  while(num_trials_ > 0) {
    mc->attempt(1);
    --num_trials_;
    DEBUG("num_trials " << num_trials_);
  }
  while(until_num_particles_ > 0 &&
        mc->configuration().num_particles() != until_num_particles_) {
    mc->attempt(1);
    DEBUG("num_particles " << mc->configuration().num_particles());
  }
  if (for_hours_ > 0) {
    const double begin = cpu_hours();
    while (for_hours_ > cpu_hours() - begin) {
      mc->attempt(10);
    }
  }
  if (until_criteria_complete_) {
    while (!mc->criteria().is_complete()) {
      mc->attempt(10);
    }
  }
}

RemoveTrial::RemoveTrial(argtype * args) {
  class_name_ = "RemoveTrial";
  index_ = integer("index", args, -1);
  name_ = str("name", args, "");
  all_ = boolean("all", args, false);
}
RemoveTrial::RemoveTrial(argtype args) : RemoveTrial(&args) {
  check_all_used(args);
}

class MapRemoveTrial {
 public:
  MapRemoveTrial() {
    auto obj = MakeRemoveTrial();
    obj->deserialize_map()["RemoveTrial"] = obj;
  }
};

static MapRemoveTrial mapper_RemoveTrial = MapRemoveTrial();

RemoveTrial::RemoveTrial(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 3854, "mismatch version: " << version);
  feasst_deserialize(&index_, istr);
  feasst_deserialize(&name_, istr);
  feasst_deserialize(&all_, istr);
}

void RemoveTrial::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(3854, ostr);
  feasst_serialize(index_, ostr);
  feasst_serialize(name_, ostr);
  feasst_serialize(all_, ostr);
}

void RemoveTrial::run(MonteCarlo * mc) {
  if (!name_.empty()) {
    for (int trial = 0; trial < mc->trials().num(); ++trial) {
      if (mc->trial(trial).class_name() == name_) {
        ASSERT(index_ < 0 || trial == index_,
          "RemoveTrial cannot specify both index and name");
        index_ = trial;
        break;
      }
    }
  }
  if (index_ >= 0) {
    mc->remove_trial(index_);
  }
  if (all_) {
    for (int i = 0; mc->trials().num() > 0; ++i) {
      mc->remove_trial(0);
      ASSERT(i < 1e5, "too many trials. Infinite loop?");
    }
  }
}

RemoveModify::RemoveModify(argtype * args) {
  class_name_ = "RemoveModify";
  index_ = integer("index", args, -1);
  name_ = str("name", args, "");
  all_ = boolean("all", args, false);
}
RemoveModify::RemoveModify(argtype args) : RemoveModify(&args) {
  check_all_used(args);
}

class MapRemoveModify {
 public:
  MapRemoveModify() {
    auto obj = MakeRemoveModify();
    obj->deserialize_map()["RemoveModify"] = obj;
  }
};

static MapRemoveModify mapper_RemoveModify = MapRemoveModify();

RemoveModify::RemoveModify(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2045, "mismatch version: " << version);
  feasst_deserialize(&index_, istr);
  feasst_deserialize(&name_, istr);
  feasst_deserialize(&all_, istr);
}

void RemoveModify::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(2045, ostr);
  feasst_serialize(index_, ostr);
  feasst_serialize(name_, ostr);
  feasst_serialize(all_, ostr);
}

void RemoveModify::run(MonteCarlo * mc) {
  DEBUG("name " << name_);
  if (!name_.empty()) {
    for (int modify = 0; modify < mc->num_modifiers(); ++modify) {
      DEBUG("mod " << mc->modify(modify).class_name());
      if (mc->modify(modify).class_name() == name_) {
        ASSERT(index_ < 0 || modify == index_,
          "RemoveModify cannot specify both index and name");
        index_ = modify;
        DEBUG("removing " << modify);
        break;
      }
    }
  }
  DEBUG("index " << index_);
  if (index_ >= 0) {
    DEBUG("removing " << index_);
    mc->remove_modify(index_);
  }
  if (all_) {
    for (int i = 0; mc->num_modifiers() > 0; ++i) {
      mc->remove_modify(0);
      ASSERT(i < 1e5, "Infinite loop?");
    }
  }
}

WriteCheckpoint::WriteCheckpoint(argtype * args) {
  class_name_ = "WriteCheckpoint";
}
WriteCheckpoint::WriteCheckpoint(argtype args) : WriteCheckpoint(&args) {
  check_all_used(args);
}

class MapWriteCheckpoint {
 public:
  MapWriteCheckpoint() {
    auto obj = MakeWriteCheckpoint();
    obj->deserialize_map()["WriteCheckpoint"] = obj;
  }
};

static MapWriteCheckpoint mapper_WriteCheckpoint = MapWriteCheckpoint();

WriteCheckpoint::WriteCheckpoint(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 5694, "mismatch version: " << version);
}

void WriteCheckpoint::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(5694, ostr);
}

void WriteCheckpoint::run(MonteCarlo * mc) {
  mc->write_checkpoint();
}

AddReference::AddReference(argtype * args) {
  class_name_ = "AddReference";
  potential_index_ = integer("potential_index", args, 0);
  cutoff_ = dble("cutoff", args, -1);
  use_cell_ = boolean("use_cell", args, false);
}
AddReference::AddReference(argtype args) : AddReference(&args) {
  check_all_used(args);
}

class MapAddReference {
 public:
  MapAddReference() {
    auto obj = MakeAddReference();
    obj->deserialize_map()["AddReference"] = obj;
  }
};

static MapAddReference mapper_AddReference = MapAddReference();

AddReference::AddReference(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 8473, "mismatch version: " << version);
  feasst_deserialize(&potential_index_, istr);
  feasst_deserialize(&cutoff_, istr);
  feasst_deserialize(&use_cell_, istr);
}

void AddReference::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(8473, ostr);
  feasst_serialize(potential_index_, ostr);
  feasst_serialize(cutoff_, ostr);
  feasst_serialize(use_cell_, ostr);
}

void AddReference::run(MonteCarlo * mc) {
  const Potential& pot = mc->system().potential(potential_index_);
  std::stringstream ss;
  pot.serialize(ss);
  std::shared_ptr<Potential> ref = std::make_shared<Potential>(ss);
  if (cutoff_ > 0) {
    const Configuration& config = mc->configuration();
    ref->set_model_params(config);
    for (int site_type = 0; site_type < config.num_site_types(); ++site_type) {
      ref->set_model_param("cutoff", site_type, cutoff_);
    }
    mc->add_to_reference(ref);
    if (use_cell_) {
      ref->set_visit_model_(MakeVisitModelCell({{"min_length",
                                                 str(cutoff_)}}));
    }
  } else {
    ASSERT(!use_cell_, "use_cell requires cutoff");
  }
}

}  // namespace feasst
