#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "monte_carlo/include/monte_carlo.h"
#include "gibbs/include/copy_following_lines.h"

namespace feasst {

FEASST_MAPPER(CopyFollowingLines,);

CopyFollowingLines::CopyFollowingLines(argtype * args) {
  class_name_ = "CopyFollowingLines";
  for_num_configurations_ = integer("for_num_configurations", args, 1);
  DEBUG("for num " << for_num_configurations_);
  std::string start;
  start.assign("replace");
  if (used(start, *args)) {
    FATAL("replace is deprecated. There seems to be no need for this usage?");
    std::vector<std::string> rpls = split(feasst::str("replace", args), ',');
    std::vector<std::string> wths = split(feasst::str("with", args), ',');
    ASSERT(rpls.size() == wths.size(), "The number of comma-separated values "
      << "given to \"replace\" must equal the number given to \"with\"");
    for (int i = 0; i < static_cast<int>(rpls.size()); ++i) {
      replace_.push_back({rpls[i], wths[i]});
    }
  } else {
    int index = 0;
    std::stringstream key;
    key << start << index;
    while (used(key.str(), *args)) {
      WARN("replace[i] is deprecated. Use comma-separated input.");
      replace_.push_back({str(key.str(), args), str("with"+str(index), args)});
      ++index;
      ASSERT(index < 1e8, "index(" << index << ") is large. Infinite loop?");
      key.str("");
      key << start << index;
    }
  }
  replace_with_index_ = str("replace_with_index", args, "");
}
CopyFollowingLines::CopyFollowingLines(argtype args) : CopyFollowingLines(&args) {
  feasst_check_all_used(args);
}

CopyFollowingLines::CopyFollowingLines(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 3406 && version <= 3407, "mismatch version: " << version);
  feasst_deserialize(&for_num_configurations_, istr);
  feasst_deserialize(&replace_, istr);
  if (version >= 3407) {
    feasst_deserialize(&replace_with_index_, istr);
  }
}

void CopyFollowingLines::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(3407, ostr);
  feasst_serialize(for_num_configurations_, ostr);
  feasst_serialize(replace_, ostr);
  feasst_serialize(replace_with_index_, ostr);
}

void CopyFollowingLines::run(MonteCarlo * mc) {
  WARN("Deprecated CopyFollowingLines. Use For instead.");
  DEBUG("setting " << for_num_configurations_);
  mc->set_parse_for_num_configs(for_num_configurations_);
  mc->set_parse_replace(replace_);
  mc->set_replace_with_index(replace_with_index_);
}

}  // namespace feasst
