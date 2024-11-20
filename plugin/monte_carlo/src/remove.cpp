#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "monte_carlo/include/analyze_factory.h"
#include "monte_carlo/include/modify_factory.h"
#include "monte_carlo/include/trial_factory.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/remove.h"

namespace feasst {

FEASST_MAPPER(Remove, argtype({{"name", "test"}}));

Remove::Remove(argtype * args) {
  class_name_ = "Remove";
  DEBUG("parse names");
  std::string start;
  start.assign("name");
  if (used(start, *args)) {
    names_.push_back(feasst::str(start, args));
  } else {
    int index = 0;
    std::stringstream key;
    key << start << index;
    while (used(key.str(), *args)) {
      names_.push_back(feasst::str(key.str(), args));
      ++index;
      ASSERT(index < 1e8, "index(" << index << ") is very high. Infinite loop?");
      key.str("");
      key << start << index;
    }
  }
  DEBUG("names: " << feasst_str(names_));

  DEBUG("parse name_contains");
  start.assign("name_contains");
  DEBUG("start " << start);
  if (used(start, *args)) {
    name_contains_.push_back(feasst::str(start, args));
  } else {
    int index = 0;
    std::stringstream key;
    key << start << index;
    while (used(key.str(), *args)) {
      name_contains_.push_back(feasst::str(key.str(), args));
      ++index;
      ASSERT(index < 1e8, "index(" << index << ") is very high. Infinite loop?");
      key.str("");
      key << start << index;
    }
  }
  DEBUG("name_contains: " << feasst_str(name_contains_));
  all_trials_ = boolean("all_trials", args, false);
  all_analyzers_ = boolean("all_analyzers", args, false);
  all_modifiers_ = boolean("all_modifiers", args, false);
  ASSERT(all_trials_ || all_analyzers_ || all_modifiers_ ||
         static_cast<int>(names_.size() + name_contains_.size()) > 0,
    "Nothing to be removed.");
}
Remove::Remove(argtype args) : Remove(&args) {
  feasst_check_all_used(args);
}

Remove::Remove(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6748, "mismatch version: " << version);
  feasst_deserialize(&names_, istr);
  feasst_deserialize(&name_contains_, istr);
  feasst_deserialize(&all_trials_, istr);
  feasst_deserialize(&all_analyzers_, istr);
  feasst_deserialize(&all_modifiers_, istr);
}

void Remove::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(6748, ostr);
  feasst_serialize(names_, ostr);
  feasst_serialize(name_contains_, ostr);
  feasst_serialize(all_trials_, ostr);
  feasst_serialize(all_analyzers_, ostr);
  feasst_serialize(all_modifiers_, ostr);
}

void Remove::run(MonteCarlo * mc) {
  DEBUG("search names");
  for (const std::string& name : names_) {
    bool removed = false;
    ASSERT(!name.empty(), "err");
    DEBUG("Search analyze names");
    for (int analyze = 0; analyze < mc->num_analyzers(); ++analyze) {
      DEBUG("an " << mc->analyze(analyze).class_name());
      if (mc->analyze(analyze).class_name() == name) {
        DEBUG("removing analyze:" << analyze);
        mc->remove_analyze(analyze);
        removed = true;
        break;
      }
    }
    if (removed) continue;
    DEBUG("Search modify names");
    for (int modify = 0; modify < mc->num_modifiers(); ++modify) {
      DEBUG("an " << mc->modify(modify).class_name());
      if (mc->modify(modify).class_name() == name) {
        DEBUG("removing modify:" << modify);
        mc->remove_modify(modify);
        removed = true;
        break;
      }
    }
    if (removed) continue;
    DEBUG("Search trial names");
    for (int trial = 0; trial < mc->trials().num(); ++trial) {
      if (mc->trial(trial).class_name() == name) {
        DEBUG("removing trial:" << trial);
        mc->remove_trial(trial);
        removed = true;
        break;
      }
    }
    if (removed) continue;
    FATAL("Nothing found to removed with the given name:" << name);
  }

  DEBUG("search name_contains");
  for (const std::string& contains : name_contains_) {
    bool removed = false;
    DEBUG("Search trial names for " << contains);
    for (int trial = mc->trials().num() - 1; trial >= 0; --trial) {
      const std::string name = mc->trial(trial).name_or_description();
      DEBUG("trial name: " << name);
      if (name.find(contains) != std::string::npos) {
        removed = true;
        mc->remove_trial(trial);
        DEBUG("removing trial: " << trial);
      }
    }
    ASSERT(removed, "Nothing found to removed which contains:" << contains);
  }
  if (all_trials_) {
    for (int trial = mc->trials().num() - 1; trial >= 0; --trial) {
      mc->remove_trial(trial);
    }
  }
  if (all_analyzers_) {
    for (int analyze = mc->analyze_factory().num() - 1; analyze >= 0; --analyze) {
      mc->remove_analyze(analyze);
    }
  }
  if (all_modifiers_) {
    for (int modify = mc->modify_factory().num() - 1; modify >= 0; --modify) {
      mc->remove_modify(modify);
    }
  }
}

}  // namespace feasst
