#include "utils/include/arguments.h"
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
  until_file_exists_ = str("until_file_exists", args, "");
  class_name_ = "Run";
}
Run::Run(argtype args) : Run(&args) {
  feasst_check_all_used(args);
}

FEASST_MAPPER(Run,);

Run::Run(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 3854 && version <= 3856, "mismatch version: " << version);
  feasst_deserialize(&num_trials_, istr);
  feasst_deserialize(&until_num_particles_, istr);
  if (version >= 3855) {
    feasst_deserialize(&configuration_index_, istr);
  }
  feasst_deserialize(&particle_type_, istr);
  feasst_deserialize(&for_hours_, istr);
  feasst_deserialize(&until_criteria_complete_, istr);
  if (version >= 3856) {
    feasst_deserialize(&until_file_exists_, istr);
  }
}

void Run::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(3856, ostr);
  feasst_serialize(num_trials_, ostr);
  feasst_serialize(until_num_particles_, ostr);
  feasst_serialize(configuration_index_, ostr);
  feasst_serialize(particle_type_, ostr);
  feasst_serialize(for_hours_, ostr);
  feasst_serialize(until_criteria_complete_, ostr);
  feasst_serialize(until_file_exists_, ostr);
}

void Run::run(MonteCarlo * mc) {
  mc->run_num_trials(num_trials_);
  mc->run_until_num_particles(until_num_particles_,
                              particle_type_,
                              configuration_index_);
  mc->run_for_hours(for_hours_);
  if (until_criteria_complete_) {
    mc->run_until_complete();
  }
  mc->run_until_file_exists(until_file_exists_);
}

}  // namespace feasst
