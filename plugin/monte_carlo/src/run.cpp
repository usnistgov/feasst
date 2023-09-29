
#include <fstream>
#include <sstream>
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
  configuration_index_ = integer("configuration_index", args, 0);
  particle_type_ = integer("particle_type", args, -1);
  for_hours_ = dble("for_hours", args, -1);
  until_criteria_complete_ = boolean("until_criteria_complete", args, false);
  class_name_ = "Run";
}
Run::Run(argtype args) : Run(&args) {
  FEASST_CHECK_ALL_USED(args);
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
  ASSERT(version >= 3854 && version <= 3855, "mismatch version: " << version);
  feasst_deserialize(&num_trials_, istr);
  feasst_deserialize(&until_num_particles_, istr);
  if (version >= 3855) {
    feasst_deserialize(&configuration_index_, istr);
  }
  feasst_deserialize(&particle_type_, istr);
  feasst_deserialize(&for_hours_, istr);
  feasst_deserialize(&until_criteria_complete_, istr);
}

void Run::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(3855, ostr);
  feasst_serialize(num_trials_, ostr);
  feasst_serialize(until_num_particles_, ostr);
  feasst_serialize(configuration_index_, ostr);
  feasst_serialize(particle_type_, ostr);
  feasst_serialize(for_hours_, ostr);
  feasst_serialize(until_criteria_complete_, ostr);
}

void Run::run(MonteCarlo * mc) {
  while (num_trials_ > 0) {
    mc->attempt(1);
    --num_trials_;
    DEBUG("num_trials " << num_trials_);
  }
  const Configuration& conf = mc->configuration(configuration_index_);
  while ((until_num_particles_ > 0) &&
         ((particle_type_ == -1 && (conf.num_particles() != until_num_particles_)) ||
          (particle_type_ != -1 && (conf.num_particles_of_type(particle_type_) != until_num_particles_)))) {
    mc->attempt(1);
    DEBUG("num_particles " << conf.num_particles());
  }
  if (for_hours_ > 0) {
    const double begin = cpu_hours();
    while (for_hours_ > cpu_hours() - begin) {
      mc->attempt(10);
    }
  }
  if (until_criteria_complete_) {
    DEBUG("num iterations " << mc->criteria().num_iterations());
    DEBUG("num iterations to complete " << mc->criteria().num_iterations_to_complete());
    DEBUG("mc->criteria().is_complete() " << mc->criteria().is_complete());
    while (!mc->criteria().is_complete()) {
      DEBUG("mc->criteria().is_complete() " << mc->criteria().is_complete());
      mc->attempt(1);
    }
  }
}

RemoveTrial::RemoveTrial(argtype * args) {
  class_name_ = "RemoveTrial";
  index_ = integer("index", args, -1);
  name_ = str("name", args, "");
  all_ = boolean("all", args, false);
  name_contains_ = str("name_contains", args, "");
}
RemoveTrial::RemoveTrial(argtype args) : RemoveTrial(&args) {
  FEASST_CHECK_ALL_USED(args);
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
  ASSERT(version >= 3854 && version <= 3855, "mismatch version: " << version);
  feasst_deserialize(&index_, istr);
  feasst_deserialize(&name_, istr);
  if (version >= 3855) {
    feasst_deserialize(&name_contains_, istr);
  }
  feasst_deserialize(&all_, istr);
}

void RemoveTrial::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(3855, ostr);
  feasst_serialize(index_, ostr);
  feasst_serialize(name_, ostr);
  feasst_serialize(name_contains_, ostr);
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
    ASSERT(index_ != -1, "No Trial of name: " << name_);
  }
  if (index_ >= 0) {
    mc->remove_trial(index_);
  }
  if (!name_contains_.empty()) {
    for (int trial = mc->trials().num() - 1; trial >= 0; --trial) {
      std::string name = mc->trial(trial).class_name();
      if (name == "Trial") {
        name = mc->trial(trial).description();
      }
      if (name.find(name_contains_) != std::string::npos) {
        mc->remove_trial(trial);
      }
    }
  }
  if (all_) {
    for (int i = 0; mc->trials().num() > 0; ++i) {
      mc->remove_trial(0);
      ASSERT(i < 1e5, "too many trials. Infinite loop?");
    }
  }
}

RemoveAnalyze::RemoveAnalyze(argtype * args) {
  class_name_ = "RemoveAnalyze";
  index_ = integer("index", args, -1);
  name_ = str("name", args, "");
  all_ = boolean("all", args, false);
}
RemoveAnalyze::RemoveAnalyze(argtype args) : RemoveAnalyze(&args) {
  FEASST_CHECK_ALL_USED(args);
}

class MapRemoveAnalyze {
 public:
  MapRemoveAnalyze() {
    auto obj = MakeRemoveAnalyze();
    obj->deserialize_map()["RemoveAnalyze"] = obj;
  }
};

static MapRemoveAnalyze mapper_RemoveAnalyze = MapRemoveAnalyze();

RemoveAnalyze::RemoveAnalyze(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 7985, "mismatch version: " << version);
  feasst_deserialize(&index_, istr);
  feasst_deserialize(&name_, istr);
  feasst_deserialize(&all_, istr);
}

void RemoveAnalyze::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(7985, ostr);
  feasst_serialize(index_, ostr);
  feasst_serialize(name_, ostr);
  feasst_serialize(all_, ostr);
}

void RemoveAnalyze::run(MonteCarlo * mc) {
  DEBUG("name " << name_);
  if (!name_.empty()) {
    for (int analyze = 0; analyze < mc->num_analyzers(); ++analyze) {
      DEBUG("an " << mc->analyze(analyze).class_name());
      if (mc->analyze(analyze).class_name() == name_) {
        ASSERT(index_ < 0 || analyze == index_,
          "RemoveAnalyze cannot specify both index and name");
        index_ = analyze;
        DEBUG("removing " << analyze);
        break;
      }
    }
    ASSERT(index_ != -1, "No Analyze of name: " << name_);
  }
  DEBUG("index " << index_);
  if (index_ >= 0) {
    DEBUG("removing " << index_);
    mc->remove_analyze(index_);
  }
  if (all_) {
    for (int i = 0; mc->num_analyzers() > 0; ++i) {
      mc->remove_analyze(0);
      ASSERT(i < 1e5, "Infinite loop?");
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
  FEASST_CHECK_ALL_USED(args);
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
    ASSERT(index_ != -1, "No Modify of name: " << name_);
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
  FEASST_CHECK_ALL_USED(args);
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

WriteStepper::WriteStepper(argtype * args) {
  class_name_ = "WriteStepper";
}
WriteStepper::WriteStepper(argtype args) : WriteStepper(&args) {
  FEASST_CHECK_ALL_USED(args);
}

class MapWriteStepper {
 public:
  MapWriteStepper() {
    auto obj = MakeWriteStepper();
    obj->deserialize_map()["WriteStepper"] = obj;
  }
};

static MapWriteStepper mapper_WriteStepper = MapWriteStepper();

WriteStepper::WriteStepper(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 7369, "mismatch version: " << version);
}

void WriteStepper::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(7369, ostr);
}

void WriteStepper::run(MonteCarlo * mc) {
  FATAL("not implemented");
}

ConvertToRefPotential::ConvertToRefPotential(argtype * args) {
  class_name_ = "ConvertToRefPotential";
  potential_index_ = integer("potential_index", args, 0);
  cutoff_ = dble("cutoff", args, -1);
  use_cell_ = boolean("use_cell", args, false);
}
ConvertToRefPotential::ConvertToRefPotential(argtype args) : ConvertToRefPotential(&args) {
  FEASST_CHECK_ALL_USED(args);
}

class MapConvertToRefPotential {
 public:
  MapConvertToRefPotential() {
    auto obj = MakeConvertToRefPotential();
    obj->deserialize_map()["ConvertToRefPotential"] = obj;
  }
};

static MapConvertToRefPotential mapper_ConvertToRefPotential = MapConvertToRefPotential();

ConvertToRefPotential::ConvertToRefPotential(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 8473, "mismatch version: " << version);
  feasst_deserialize(&potential_index_, istr);
  feasst_deserialize(&cutoff_, istr);
  feasst_deserialize(&use_cell_, istr);
}

void ConvertToRefPotential::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(8473, ostr);
  feasst_serialize(potential_index_, ostr);
  feasst_serialize(cutoff_, ostr);
  feasst_serialize(use_cell_, ostr);
}

void ConvertToRefPotential::run(MonteCarlo * mc) {
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

RefPotential::RefPotential(argtype * args) {
  class_name_ = "RefPotential";
  reference_index_ = integer("reference_index", args, 0);
  args_ = *args;
  args->clear();
}
RefPotential::RefPotential(argtype args) : RefPotential(&args) {
  FEASST_CHECK_ALL_USED(args);
}

class MapRefPotential {
 public:
  MapRefPotential() {
    auto obj = MakeRefPotential();
    obj->deserialize_map()["RefPotential"] = obj;
  }
};

static MapRefPotential mapper_RefPotential = MapRefPotential();

RefPotential::RefPotential(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 4017, "mismatch version: " << version);
  feasst_deserialize(&reference_index_, istr);
  feasst_deserialize(&args_, istr);
}

void RefPotential::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(4017, ostr);
  feasst_serialize(reference_index_, ostr);
  feasst_serialize(args_, ostr);
}

void RefPotential::run(MonteCarlo * mc) {
  mc->add_to_reference(MakePotential(args_), reference_index_);
}

WriteModelParams::WriteModelParams(argtype * args) {
  class_name_ = "WriteModelParams";
  file_name_ = str("file_name", args);
  potential_index_ = integer("potential_index", args, -1);
  reference_index_ = integer("reference_index", args, -1);
}
WriteModelParams::WriteModelParams(argtype args) : WriteModelParams(&args) {
  FEASST_CHECK_ALL_USED(args);
}

class MapWriteModelParams {
 public:
  MapWriteModelParams() {
    auto obj = MakeWriteModelParams({{"file_name", "place_holder"}});
    obj->deserialize_map()["WriteModelParams"] = obj;
  }
};

static MapWriteModelParams mapper_WriteModelParams = MapWriteModelParams();

WriteModelParams::WriteModelParams(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2890, "mismatch version: " << version);
  feasst_deserialize(&file_name_, istr);
  feasst_deserialize(&potential_index_, istr);
  feasst_deserialize(&reference_index_, istr);
}

void WriteModelParams::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(2890, ostr);
  feasst_serialize(file_name_, ostr);
  feasst_serialize(potential_index_, ostr);
  feasst_serialize(reference_index_, ostr);
}

void WriteModelParams::run(MonteCarlo * mc) {
  std::ofstream file(file_name_);
  const Configuration& config = mc->configuration();
  if (potential_index_ == -1) {
    file << config.model_params().str();
  } else {
    if (reference_index_ == -1) {
      const Potential& potential = mc->system().potential(potential_index_);
      file << potential.model_params(config).str();
    } else {
      const Potential& potential = mc->system().reference(reference_index_,
                                                          potential_index_);
      file << potential.model_params(config).str();
    }
  }
}

}  // namespace feasst
