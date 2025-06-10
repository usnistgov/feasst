#include <iostream>
#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "monte_carlo/include/monte_carlo.h"
#include "gibbs/include/copy_next_line.h"

namespace feasst {

FEASST_MAPPER(CopyNextLine,);

CopyNextLine::CopyNextLine(argtype * args) {
  class_name_ = "CopyNextLine";
  std::string start;
  start.assign("replace");
  if (used(start, *args)) {
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
}
CopyNextLine::CopyNextLine(argtype args) : CopyNextLine(&args) {
  feasst_check_all_used(args);
}

CopyNextLine::CopyNextLine(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 5345 && version <= 5345, "mismatch version: " << version);
  feasst_deserialize(&replace_, istr);
}

void CopyNextLine::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(5345, ostr);
  feasst_serialize(replace_, ostr);
}

void CopyNextLine::run(MonteCarlo * mc) {
  WARN("Deprecated CopyNextLine. Use For instead.");
  std::pair<std::string, argtype> next_arg = mc->next_arg();
  for (const std::vector<std::string>& rep : replace_) {
    str(rep[0], &(next_arg.second)); // remove arg
    next_arg.second.insert({rep[0], rep[1]});
  }
  arglist newline = {next_arg};
  std::cout << "# Copied line: " << next_arg.first << " " << str(next_arg.second) << std::endl;
  mc->parse_args(&newline, true);  // parse silently so not in reproducible logs
}

}  // namespace feasst
